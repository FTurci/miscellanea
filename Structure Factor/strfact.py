import numpy as np
from numpy.fft import fftn, fftshift
import pylab as pl

def radial_profile(data, center):
    y, x,z = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2+ (z - center[2])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

# print i
conf=np.loadtxt("conf.xyz", skiprows=2, usecols=[1,2,3])
Natoms=conf.shape[0]
Ndiv=170

# x = np.linspace(-L, L, Ndiv+1, endpoint=True)
f, edges = np.histogramdd(conf, bins=(Ndiv, Ndiv, Ndiv))

ftk = (fftshift(fftn(fftshift(f))))
sk = np.abs(ftk**2) / float(Natoms)
# find max
centre=np.unravel_index(sk.argmax(),sk.shape)
print centre
p=radial_profile(sk, centre)
pl.loglog(p)
pl.show()
