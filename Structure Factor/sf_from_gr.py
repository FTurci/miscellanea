import numpy as np

# path to g_r file
filename="data2.txt"
# read r,g
L=63.9719
N=10000
rho=N/L**3
# bin,r,g=np.loadtxt(filename, unpack=True)
r,g=np.loadtxt(filename, unpack=True)

dq=0.005
qs=np.arange(2*np.pi/r.max(),100,dq)
# computing the Structure factor as an isotropic Fourier Transform of the g(r)
S=1.+4.*np.pi*rho*np.array([np.trapz((g-1.)*r*np.sin(q*r)/q,x=r) for q in qs])
# plotting
import pylab as pl
pl.clf()
pl.plot(qs,S)
pl.show()