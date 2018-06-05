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
from functions import getgamma,getkappa,getpowerspectrum,getpowerspectra,getfields,plotting,plotpowerspectra,getlineprofiles,getxsquare

pointnum = 100
fieldnum = 10
fieldsize = 1. #fieldsize in degree
mapsave = "n"
lineplot = 'n'
numruns = 0



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
#plt.show()
pskappa,pskappa2,pskappa2r,pskappa2i,pskappadiff = getpowerspectra(kappa,kappa2back)
function,function2,function3 = getlineprofiles(kappa,kappa2back,fieldnum,pointnum,weightf,False,20)
for i in range(numruns):
    kappa,gamma,kappaback,kappa2back,weightf = getfields(gesnum,fieldnum,pointnum,Dhat,Dstarhat)
    pskappat,pskappa2t,pskappa2rt,pskappa2it,pskappadifft = getpowerspectra(kappa,kappa2back)
    functiont,function2t,function3t = getlineprofiles(kappa,kappa2back,fieldnum,pointnum,weightf,False,20)
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

function = function/(numruns+1)
function2 = function2/(numruns+1)
function3 = function3/(numruns+1)


ps1 = getpowerspectrum(kappa)
ps2 = getpowerspectrum(kappa2back)
ps3 = getpowerspectrum(gamma)
ps4 = getpowerspectrum(gamma2)
ps5 = getpowerspectrum(kappa2back-kappa)
ps6 = getpowerspectrum(kappa2back.real - kappa)

def getfouriers(kappa,kappa2back,gamma,weightf,fieldnum,gesnum,pointnum):
    gamma2 = makefield.applydepthc(gamma,weightf,fieldnum,gesnum,pointnum)
    pp1 = np.abs(np.fft.fftshift(np.fft.fft2(kappa)))**2
    pp2 = np.abs(np.fft.fftshift(np.fft.fft2(kappa2back.real-kappa)))**2
    pp3 = np.abs(np.fft.fftshift(np.fft.fft2(kappa2back.imag)))**2
    pp4 = np.abs(np.fft.fftshift(np.fft.fft2(gamma)))**2
    pp5 = np.abs(np.fft.fftshift(np.fft.fft2(gamma2)))**2

    return pp1,pp2,pp3,pp4,pp5

def plotfouriers(pp1,pp2,pp3,pp4,pp5,gesnum,fieldnum):
    xsquare = getxsquare(pp1)
    ticknum = 10
    tickpositions = np.arange(ticknum+1)*gesnum/ticknum
#    plt.figure(1)
    plt.imshow(np.log10(pp1))
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$\\log10(P_{\\kappa})$')
    plt.savefig('pp1_1.png')
    plt.close()

#    plt.figure(2)
    plt.imshow(np.log10(pp2))
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$\\log10(P_{\\kappa2_r-\\kappa})$')
    plt.savefig('pp2_1.png')
    plt.close()


#    plt.figure(3)
    plt.imshow(np.log10(pp3))
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$\\log10(P_{\\kappa2_i})$')
    plt.savefig('pp3_1.png')
    plt.close()

#    plt.figure(4)
    plt.imshow(np.log10(pp4))
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$\\log10(P_{\\gamma})$')
    plt.savefig('pp4_1.png')
    plt.close()

#    plt.figure(5)
    plt.imshow(np.log10(pp5))
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$\\log10(P_{\\gamma2})$')
    plt.savefig('pp5_1.png')
    plt.close()

#    plt.figure(1)
    plt.imshow(xsquare*pp1)
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$l^2P_{\\kappa}$')
    plt.savefig('pp1_2.png')
    plt.close()

#    plt.figure(2)
    plt.imshow(xsquare*pp2)
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$l^2P_{\\kappa2_r-\\kappa}$')
    plt.savefig('pp2_2.png')
    plt.close()

#    plt.figure(3)
    plt.imshow(xsquare*pp3)
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$l^2P_{\\kappa2_i}$')
    plt.savefig('pp3_2.png')
    plt.close()

#    plt.figure(4)
    plt.imshow(xsquare*pp4)
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$l^2P_{\\gamma}$')
    plt.savefig('pp4_2.png')
    plt.close()

#    plt.figure(5)
    plt.imshow(xsquare*pp5)
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$l^2P_{\\gamma2}$')
    plt.savefig('pp5_2.png')
    plt.close()


x = np.array(range(len(ps1)))
plt.figure(1)
plt.plot(x,x*x*ps1,label = 'kappa')
plt.plot(x,x*x*ps3,label = 'gamma')
plt.plot(x,x*x*ps2,label = 'kappa2')
plt.plot(x,x*x*ps4,label = 'gamma2')
plt.legend()

plt.figure(2)
plt.plot(x,x*x*ps5,label = '1')
plt.plot(x,x*x*ps6,label = '2')
plt.legend()

plt.show()

y = np.array(range(len(pskappa)))
z = y/fieldnum

plt.plot(z,z*z*fieldnum**2*pskappa2i)
plt.show()
