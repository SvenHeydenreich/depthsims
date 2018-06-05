import numpy as np
from scipy import fftpack
import radialProfile

sigma8 = 0.02
def depthfunction():
    return np.random.normal(1,0.1)


def applydepth(a,fieldnum,num,pointnum):
    b = np.zeros((num,num),dtype='complex')
    u = np.zeros((fieldnum,fieldnum))
    for i in range(fieldnum):
        for j in range(fieldnum):
                u[i,j] = depthfunction()
    for i in range(num):
        for j in range(num):
            b[i,j]=a[i,j]*u[int(i/pointnum),int(j/pointnum)]
    return b,u

def makefield(num,fieldnum,pointnum):
    print('Generating random field')
    ghat = np.zeros((num,num))
    num1 = num-1
    for i in range(num):
        for j in range(num):
            ghat[i,j] = np.random.normal(0,1/(sigma8*np.sqrt((i-num1/2)**2+(j-num1/2)**2)))



    ghat2 = fftpack.ifftshift(ghat)
    g = fftpack.ifft2(ghat2)
    g = g.real
    g = 20.*g + 1

    return g

def getpowerspectrum(g):
    gback = fftpack.fft2(g)
    gback2 = fftpack.fftshift(gback)
    gback2 = np.abs(gback2)**2
    ps2 = radialProfile.azimuthalAverage(gback2)
    return ps2

#plt.figure(1)
#plt.imshow(g)
#plt.figure(2)
#plt.imshow(g2)
