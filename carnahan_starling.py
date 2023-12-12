#!/usr/local/bin/python
import sys
import numpy
from numpy import pi
import numpy as np
from scipy.stats import moment
# from scipy.integrate import 
 
def discrete_std(a,b):
  return np.sqrt(((b-a+1)**2-1.)/12.)

def cs(eta):
	return (1+eta+eta**2-eta**3)/(1-eta)**3*eta*6./pi

def mcs(eta,polydispersity ):
	n = 5
	p = 1.0/n
	mean = 1.0
	diams = mean+np.arange(-(n-1)/2, (n+1)/2)/discrete_std(0,n-1)*polydispersity/mean
	radii = diams/2.0

	population = np.random.choice(radii, size=100000)

	m1 = mean
	m2 = moment(population,2)+m1**2
	m3 = moment(population,3)+m1**3

	print(m1,m2,m3)
	print ("var vs moment",population.var(), moment(population,2) )

	O1 = m1*m2/m3
	O2 = m2**3/m3**2
	return 1/(1-eta)+O1*3*eta/(1-eta)**2+O2*eta**2*(3-eta)/(1-eta)**3	

if __name__ == "__main__":
	nopt=len(sys.argv)
	if nopt<2:
		print ("\n!! Provide an input packing fraction !!")
	if nopt==2:
	    P=cs(float(sys.argv[1]))

	    print ("CS",P)
	    print("MCS", mcs(float(sys.argv[1]), 0.08 ))
	if nopt==3:
		import pylab as pl
		etas=pl.linspace(float(sys.argv[1]),float(sys.argv[2]),100)
		pl.plot(etas, cs(etas))
		pl.show()

# cs(0.57)

