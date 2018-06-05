from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy import signal,ndimage
#import fields
import makefield
import distributegalaxies
import shearcorrelation
from astropy.convolution import convolve,convolve_fft
from time import time
import os
import subprocess
import shutil

dmax = 1500
dsteps = 750
pointnum = 400
fieldnum = 10
fieldsize = 30. #fieldsize in arcseconds
smoothingscale = 2
mapsave = "n"
galaxynum = 20000



gesnum = pointnum*fieldnum
kernelnum1 = int(gesnum)
kernelnum = kernelnum1-1
fieldsize = fieldsize/3600.
startt = time()
print("Filling the D-functions")
Dstar1,D1 = makefield.filld(kernelnum,fieldsize,pointnum)
Dstar1 = Dstar1/np.sum(np.abs(Dstar1))
D1 = D1/np.sum(np.abs(D1))
Dstar1r = Dstar1.real
Dstar1i = Dstar1.imag
D1r = D1.real
D1i = D1.imag
print("Finished filling. Time:",round(time()-startt),"s")

#print (np.sum(Dstar1))
#print (np.sum(D1))
def getgamma(kappa):
    kappa = kappa.astype("complex")
    print("Convolving kappa to gamma...")
    kappar = kappa.real
    kappai = kappa.imag
    gammar = convolve_fft(kappar,D1r,boundary="wrap") - convolve_fft(kappai,D1i,boundary="wrap")
    gammai = convolve_fft(kappar,D1i,boundary="wrap") + convolve_fft(kappai,D1r,boundary="wrap")
#    gammages = gammar + gammai*1j
    print("Smoothing...")
    gammar2 = ndimage.gaussian_filter(gammar,smoothingscale,mode="wrap")
    gammai2 = ndimage.gaussian_filter(gammai,smoothingscale,mode="wrap")
    gamma2 = gammar2 + gammai2*1.j
    return gamma2

def getkappa(gamma):
    gamma = gamma.astype("complex")
    print("Convolving gamma to kappa...")
    gammar = gamma.real
    gammai = gamma.imag
    kappar = convolve_fft(gammar,Dstar1r,boundary="wrap") - convolve_fft(gammai,Dstar1i,boundary="wrap")
    kappai = convolve_fft(gammar,Dstar1i,boundary="wrap") + convolve_fft(gammai,Dstar1r,boundary="wrap")
#    gammages = gammar + gammai*1j
    print("Smoothing...")
    kappar2 = ndimage.gaussian_filter(kappar,smoothingscale,mode="wrap")
    kappai2 = ndimage.gaussian_filter(kappai,smoothingscale,mode="wrap")
    kappa2 = kappar2 + kappai2*1.j
    return kappa2


def analyze(counter):
    print('Starting run ',counter)
    kappa = makefield.create(gesnum,np.random.randint(64000))

    gamma = getgamma(kappa)
    u,w = makefield.getdepth(fieldnum)
    xposs1,yposs1,shears1 = distributegalaxies.distributefullweighted(u,w,0.3, gamma, gesnum, pointnum, fieldnum, galaxynum, np.random.randint(64000))
    xposs2,yposs2,shears2 = distributegalaxies.distributeuniform(gamma,0.3,gesnum,galaxynum,np.random.randint(64000))

    #startt1 = time()
    print('Starting shearcorr1')
    shears1 = np.array(shears1)
    shears1 = 1000.*shears1
    shears2 = np.array(shears2)
    shears2 = 1000.*shears2

    #xip1,xim1,xix1 = shearcorrelation.analyze(xposs1,yposs1,shears1,dmax,dsteps)
    #print(int(time()-startt1))
    #print('Starting shearcorr2')
#    xip2,xim2,xix2 = shearcorrelation.analyze(xposs2,yposs2,shears2,dmax,dsteps)
    #print('finished shearcorr')

    startt2 = time()
    fil = open('galaxies.txt','w')
    for i in range(galaxynum):
        fil.write(str(xposs1[i])+' '+str(yposs1[i])+' '+str(shears1[i].real)+' '+str(shears1[i].imag)+'\n')
    #    print(shears1[i].imag)
    fil.close()
    subprocess.call("./a.out")
    print('Time: ',int(time()-startt2))
    shutil.copy('functions.txt','functions/'+str(counter)+'.txt')
    subprocess.call(['rm','functions.txt'])
    subprocess.call(['rm','galaxies.txt'])

    startt2 = time()
    fil = open('galaxies.txt','w')
    for i in range(galaxynum):
        fil.write(str(xposs2[i])+' '+str(yposs2[i])+' '+str(shears2[i].real)+' '+str(shears2[i].imag)+'\n')
    #    print(shears1[i].imag)
    fil.close()
    subprocess.call("./a.out")
    print('Time: ',int(time()-startt2))
    shutil.copy('functions.txt','functions/'+str(counter)+'uni.txt')
    subprocess.call(['rm','functions.txt'])
    subprocess.call(['rm','galaxies.txt'])


    return 0

if not(os.path.exists('functions')):
    os.makedirs('functions')

for j in range(100):
    analyze(j)
