from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy import signal,ndimage,fftpack
#import fields
import makefield
import distributegalaxies
import shearcorrelation
from astropy.convolution import convolve,convolve_fft
from time import time
import os
import subprocess
import shutil
import radialProfile
from functions import getgamma,getkappa,getpowerspectrum,getpowerspectra,getfields,plotting,plotpowerspectra,getlineprofiles,getfouriers,plotfouriers

starttime = time()
pointnum = 101
fieldnum = 20
fieldsize = 1. #fieldsize in degree
mapsave = "n"
lineplot = 'n'
numruns = 10



gesnum = pointnum*fieldnum

kernelnum1 = int(gesnum)
kernelnum = kernelnum1
fieldsize = fieldsize*np.pi/180 #fieldsize in rad
pixsize = fieldsize/pointnum #pixelsize in rad
invpixsize = 1./pixsize
print(pixsize,invpixsize)
startt = time()
print("Filling the D-functions")
Dhat,Dstarhat = makefield.filldhat(kernelnum,invpixsize,pointnum)
print("Finished filling. Time:",round(time()-startt),"s")

#print (np.sum(Dstar1))
#print (np.sum(D1))

kappa,gamma,kappaback,kappa2back,weightf = getfields(gesnum,fieldnum,pointnum,Dhat,Dstarhat)
pskappa,pskappa2,pskappa2r,pskappa2i,pskappadiff = getpowerspectra(kappa,kappa2back)
function,function2,function3 = getlineprofiles(kappa,kappa2back,fieldnum,pointnum,weightf,False,20)
pp1,pp2,pp3,pp4,pp5,pp6 = getfouriers(kappa,kappa2back,gamma,weightf,fieldnum,gesnum,pointnum)
for i in range(numruns):
    kappa,gamma,kappaback,kappa2back,weightf = getfields(gesnum,fieldnum,pointnum,Dhat,Dstarhat)
    pskappat,pskappa2t,pskappa2rt,pskappa2it,pskappadifft = getpowerspectra(kappa,kappa2back)
    functiont,function2t,function3t = getlineprofiles(kappa,kappa2back,fieldnum,pointnum,weightf,False,20)
    pp1t,pp2t,pp3t,pp4t,pp5t,pp6t = getfouriers(kappa,kappa2back,gamma,weightf,fieldnum,gesnum,pointnum)

    pp1 = pp1+pp1t
    pp2 = pp2+pp2t
    pp3 = pp3+pp3t
    pp4 = pp4+pp4t
    pp5 = pp5+pp5t
    pp6 = pp6+pp6t

    function = function + functiont
    function2 = function2 + function2t
    function3 = function3 + function3t

    pskappadiff = pskappadiff+pskappadifft
    pskappa = pskappa+pskappat
    pskappa2 = pskappa2+pskappa2t
    pskappa2r = pskappa2r+pskappa2rt
    pskappa2i = pskappa2i+pskappa2it

pskappa = pskappa/(numruns+1)
pskappa2 = pskappa2/(numruns+1)
pskappa2r = pskappa2r/(numruns+1)
pskappa2i = pskappa2i/(numruns+1)
pskappadiff = pskappadiff/(numruns+1)

pp1 = pp1/(numruns+1)
pp2 = pp2/(numruns+1)
pp3 = pp3/(numruns+1)
pp4 = pp4/(numruns+1)
pp5 = pp5/(numruns+1)
pp6 = pp6/(numruns+1)

function = function/(numruns+1)
function2 = function2/(numruns+1)
function3 = function3/(numruns+1)

plotpowerspectra(pskappa,pskappa2r,pskappa2i,pskappadiff,fieldnum)
plotting(kappa,gamma,kappaback,kappa2back,weightf,"y","test.png",pointnum,gesnum,fieldnum)
plotting(kappa,gamma,kappaback,kappa2back,weightf,"n","test2.png",pointnum,gesnum,fieldnum)
plotfouriers(pp1,pp2,pp3,pp4,pp5,pp6,gesnum,fieldnum)



mean = np.sqrt(np.mean(kappa**2))
x = range(len(function))
x = np.array(x)*pixsize*180/np.pi*60
plt.plot(x,np.sqrt(function)/mean,label = 'B')
plt.plot(x,np.sqrt(function2)/mean,label = 'E')
plt.plot(x,np.sqrt(function3)/mean,label = 'EW')
plt.legend()
#plt.yscale('log')
plt.ylim(0)
plt.xlabel('$\\theta$ [arcmin]')
plt.ylabel('$\\sigma_{E/B}/\\sigma_{\\kappa}$')
plt.savefig('lineprofiles.png')
plt.close()



np.save('results',[pskappa,pskappa2,pskappa2r,pskappa2i,pskappadiff,function,function2,function3])
print('Duration: ',np.round(time()-starttime))
