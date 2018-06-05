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

def getgamma(kappa,Dhat):
    kappa = kappa.astype("complex")
    print("Convolving kappa to gamma...")
    kappahat = np.fft.fft2(kappa)
    gammahat = kappahat * Dhat
    gamma = np.fft.ifft2(gammahat)
#    print("Smoothing...")
#    gammar2 = ndimage.gaussian_filter(gammar,smoothingscale,mode="wrap")
#    gammai2 = ndima100ge.gaussian_filter(gammai,smoothingscale,mode="wrap")
#    gamma2 = gammar2 + gammai2*1.j
    return gamma/np.pi

def getkappa(gamma,Dstarhat):
    gamma = gamma.astype("complex")
    print("Convolving gamma to kappa...")
    gammahat = np.fft.fft2(gamma)
    kappahat = gammahat*Dstarhat
    kappa = np.fft.ifft2(kappahat)
#    print("Smoothing...")
#    kappar2 = ndimage.gaussian_filter(kappar,smoothingscale,mode="wrap")
#    kappai2 = ndimage.gaussian_filter(kappai,smoothingscale,mode="wrap")
#    kappa2 = kappar2 + kappai2*1.j
    return kappa/np.pi

def getpowerspectrum(g):
    ghat = np.fft.fft2(g)
    ghatshift = fftpack.fftshift(ghat)
#    ghatshift = ghat
    gabs = np.abs(ghatshift)**2
    ps2 = radialProfile.azimuthalAverage(gabs)
    return ps2

def getfields(gesnum,fieldnum,pointnum,Dhat,Dstarhat):
    kappa = makefield.create(gesnum,np.random.randint(64000))
    gamma = getgamma(kappa,Dhat)
    weightf = makefield.drawdepth(fieldnum,gesnum,pointnum)
    gamma2 = makefield.applydepthc(gamma,weightf,fieldnum,gesnum,pointnum)
    kappaback = getkappa(gamma,Dstarhat)
    kappa2back = getkappa(gamma2,Dstarhat)
    print('Weighted:')
    print(np.sum(np.abs(kappa - kappa2back))/np.sum(np.abs(kappa)))
    print(np.sum(np.abs(kappa - kappa2back.real))/np.sum(np.abs(kappa)))
    print(np.sum(np.abs(kappa2back.imag))/np.sum(np.abs(kappa)))

    print('Uniform:')
    print(np.sum(np.abs(kappa - kappaback))/np.sum(np.abs(kappa)))
    print(np.sum(np.abs(kappa - kappaback.real))/np.sum(np.abs(kappa)))
    print(np.sum(np.abs(kappaback.imag))/np.sum(np.abs(kappa)))


    return kappa,gamma,kappaback,kappa2back,weightf

def getpowerspectra(kappa,kappa2back):
    pskappa = getpowerspectrum(kappa)
    pskappa2 = getpowerspectrum(kappa2back)
    pskappa2r = getpowerspectrum(kappa2back.real)
    pskappa2i = getpowerspectrum(kappa2back.imag)
    pskappadiff = getpowerspectrum((kappa2back.real-kappa))
    return pskappa,pskappa2,pskappa2r,pskappa2i,pskappadiff

def plotting(kappa,gamma,kappaback,kappa2back,weightf,lineplot,savename,pointnum,gesnum,fieldnum):
    lines = np.linspace(pointnum,gesnum-pointnum,fieldnum-1)


    kappamin = np.min([np.min(kappa2back.real),np.min(kappa)])
    kappamax = np.max([np.max(kappa2back.real),np.max(kappa)])

    divmin = np.min([np.min(kappa2back.real-kappa),np.min(kappa2back.imag)])
    divmax = np.max([np.max(kappa2back.real-kappa),np.max(kappa2back.imag)])

    plt.figure(1,figsize=(14, 12), dpi=300, facecolor='w', edgecolor='k')

    plt.subplot(321)
    plt.imshow(10**4*(kappa),vmin=10**4*kappamin,vmax=10**4*kappamax)
    plt.title('Kappa')
    plt.axis('off')
    plt.colorbar()
    if(lineplot=='y'):
        plt.hlines(lines,0,gesnum,linestyles='dotted')
        plt.vlines(lines,0,gesnum,linestyles='dotted')


    plt.subplot(322)
    plt.imshow(10**4*(kappa2back.real),vmin=10**4*kappamin,vmax=10**4*kappamax)
    plt.title('Re(Kappa2)')
    plt.axis('off')
    plt.colorbar()
    if(lineplot=='y'):
        plt.hlines(lines,0,gesnum,linestyles='dotted')
        plt.vlines(lines,0,gesnum,linestyles='dotted')



    plt.subplot(323)
    plt.imshow(10**4*(makefield.applydepth(kappa2back.real,1/weightf,fieldnum,gesnum,pointnum) - kappa),vmin=10**4*divmin,vmax=10**4*divmax)
    plt.title('Re(Kappa2)[1/weighted] - Kappa')
    plt.axis('off')
    plt.colorbar()
    if(lineplot=='y'):
        plt.hlines(lines,0,gesnum,linestyles='dotted')
        plt.vlines(lines,0,gesnum,linestyles='dotted')



    plt.subplot(324)
    plt.imshow(10**4*(kappa2back.imag),vmin=10**4*divmin,vmax=10**4*divmax)
    plt.title('Im(Kappa2)')
    plt.axis('off')
    plt.colorbar()
    if(lineplot=='y'):
        plt.hlines(lines,0,gesnum,linestyles='dotted')
        plt.vlines(lines,0,gesnum,linestyles='dotted')


    plt.subplot(325)
    plt.imshow(10**4*(kappaback.real - kappa))
    plt.title('Numerical Errors')
    plt.axis('off')
    plt.colorbar()

    plt.subplot(326)
    plt.imshow(weightf,interpolation='None',cmap='Greys')
    plt.title('Weighting Function')
    plt.axis('off')
    plt.colorbar()




    #plt.gca().yaxis.set_minor_formatter(NullFormatter())
    #plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,wspace=0.0)
    plt.savefig(savename)
    plt.close()

def plotpowerspectra(pskappa,pskappa2r,pskappa2i,pskappadiff,fieldnum):
    plt.close()
    x = np.array(range(len(pskappa)))/fieldnum
#    plt.figure(2,figsize=(14, 12), dpi=300, facecolor='w', edgecolor='k')
#    plt.subplot(221)
    plt.plot(x,x*x*fieldnum**2*pskappa)
    plt.ylabel('$l^2P_{\\kappa}(l)$')
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.savefig('pkappa.png')
    plt.close()
#    plt.axis('off')

#    plt.subplot(222)
    plt.plot(x,x*x*fieldnum**2*(pskappadiff))
    plt.ylabel('$l^2P_{\\kappa 2_r-\\kappa}(l)$')
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.savefig('pkappadiff.png')
    plt.close()
#    plt.axis('off')

#    plt.subplot(223)
    plt.plot(x,x*x*fieldnum**2*pskappa2r)
    plt.ylabel('$l^2P_{\\kappa 2_r}(l)$')
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.savefig('pkappa2r.png')
    plt.close()
#    plt.axis('off')

#    plt.subplot(224)
    plt.plot(x,x*x*fieldnum**2*pskappa2i)
    plt.ylabel('$l^2P_{\\kappa 2_i}(l)$')
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.savefig('pkappa2i.png')
    plt.close()
#    plt.axis('off')

#    plt.savefig('powerspectra.png')
#    plt.close()

def getlineprofile(field,fieldnum,pointnum,excludeboundary,boundary = 0):
    gesnum = fieldnum*pointnum
    function = getoneline(field,fieldnum,pointnum,excludeboundary,boundary)
    for i in range(3):
        function = function + getoneline(np.rot90(field,i+1),fieldnum,pointnum,excludeboundary,boundary)
    return function



def getoneline(field,fieldnum,pointnum,excludeboundary,boundary):
    gesnum = fieldnum*pointnum
    lines = np.linspace(0,gesnum,fieldnum+1)
    function = np.zeros(int(pointnum/2))
    if(excludeboundary):
        for i in range(len(lines)-1):
            a = field[int(lines[i]):int((lines[i+1]-lines[i])/2),:]
            b = np.sum(a,1)/gesnum
            function = function+b
    else:
        for i in range(len(lines)-1):
            for j in range(len(lines)-1):
                a = field[int(lines[i]):int(lines[i])+int((lines[i+1]-lines[i])/2),int(lines[j])+boundary:int(lines[j+1])-boundary]
                b = np.sum(a,1)
                b = b/(pointnum-2*boundary)
#                print(b.shape)
                function = function+b
        function = function/fieldnum
    return function

def getlineprofiles(kappa,kappa2back,fieldnum,pointnum,weightf,excludeboundary,boundary=20):
    gesnum = fieldnum*pointnum
    function = getlineprofile((kappa2back.imag)**2,fieldnum,pointnum,excludeboundary,boundary)
    function2 = getlineprofile((kappa2back.real-kappa)**2,fieldnum,pointnum,excludeboundary,boundary)
    field = makefield.applydepth(kappa2back.real,1/weightf,fieldnum,gesnum,pointnum) - kappa
    function3 = getlineprofile((field)**2,fieldnum,pointnum,excludeboundary,boundary)
    return function,function2,function3


def getxsquare(pp1):
    y,x = np.indices(pp1.shape)
    center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    return r**2

def getfouriers(kappa,kappa2back,gamma,weightf,fieldnum,gesnum,pointnum):
    gamma2 = makefield.applydepthc(gamma,weightf,fieldnum,gesnum,pointnum)
    kappaw = makefield.applydepth(kappa2back.real,1/weightf,fieldnum,gesnum,pointnum) - kappa
    pp1 = np.abs(np.fft.fftshift(np.fft.fft2(kappa)))**2
    pp2 = np.abs(np.fft.fftshift(np.fft.fft2(kappa2back.real-kappa)))**2
    pp3 = np.abs(np.fft.fftshift(np.fft.fft2(kappa2back.imag)))**2
    pp4 = np.abs(np.fft.fftshift(np.fft.fft2(gamma)))**2
    pp5 = np.abs(np.fft.fftshift(np.fft.fft2(gamma2)))**2
    pp6 = np.abs(np.fft.fftshift(np.fft.fft2(kappaw)))**2

    return pp1,pp2,pp3,pp4,pp5,pp6

def plotfouriers(pp1,pp2,pp3,pp4,pp5,pp6,gesnum,fieldnum):
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

    plt.imshow(np.log10(pp6))
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$\\log10(P_{\\kappa2_{r,w}-\\kappa})$')
    plt.savefig('pp6_1.png')
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

    plt.imshow(xsquare*pp6)
    plt.xticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.yticks(tickpositions,((tickpositions - gesnum/2)/fieldnum))
    plt.xlabel('$l$ [1/fieldnumber]')
    plt.ylabel('$l$ [1/fieldnumber]')
    plt.colorbar()
    plt.title('$l^2P_{\\kappa2_{r,w}-\\kappa}$')
    plt.savefig('pp6_2.png')
    plt.close()
