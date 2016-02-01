import numpy as np
import pylab as pl

L=1
N=100
D=2
random_points=np.random.uniform(0,L,size=(N,D)) # N points in D dimensions

pl.scatter(random_points[:,0], random_points[:,1])

H, edges=np.histogramdd(random_points, bins=(np.arange(0,L,L/10.),np.arange(0,L,L/10.)))

pl.figure()

# count the number of covered bins
print np.sum(H>0)
pl.show()