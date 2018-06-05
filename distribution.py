import numpy as np
import matplotlib.pyplot as plt

def getvals(zmin):
    folder = 'distributions'
    name1 = 'KiDS_2017-05-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB'
    name2 = '_rlim_percentile'
    name3 = '_blindB_Nz.asc'
    zmin = np.round(zmin,1)
    zmax = zmin+0.2
    if(zmin==0.9):
        zmax = 1.2
    inter = str(zmin).replace('.','p')+'t'+str(zmax).replace('.','p')

    vals = []
    for i in range(10):
        fil = open(folder+'/'+name1+inter+name2+str(10*i)+'t'+str(10*(i+1))+name3)
        values = fil.readlines()
        fil.close()
        l = len(values)
        zs = np.zeros(l)
        vs = np.zeros(l)
        for j in range(l):
            zs[j],vs[j]=values[j].split()
        if(i==1):
            redshifts = zs
        vals.append(vs)
    return redshifts,vals

def getredshift(values,redshifts,length):
    valmax = values.sum()
    valinter = np.zeros(length+1)
    valtemp = 0
    for i in range(length):
        valtemp = valtemp+values[i]
        valinter[i+1] = valtemp
    r = np.random.uniform(0,valmax)
    for i in range(length):
        if(valinter[i]<r and r<=valinter[i+1]):
            binnum = i
    if(binnum != length-1):
        return np.random.uniform(redshifts[binnum],redshifts[binnum+1])
    else:
        return redshifts[binnum]
