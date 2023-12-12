import pylab as pl
from numpy import pi,vectorize, zeros,column_stack
from numba import autojit,jit



@jit(nopython=True)
def spherical_FT(f,k,r,dr):
	ft=pl.zeros(len(k))
	for i in range(len(k)):
		ft[i]=4.*pi*pl.sum(r*pl.sin(k[i]*r)*f*dr)/k[i]
	return ft
@jit(nopython=True)
def inverse_spherical_FT(ff,k,r,dk):
	inft=pl.zeros(len(r))
	for i in range(len(r)):
		inft[i]=pl.sum(k*pl.sin(k*r[i])*ff*dk)/r[i]/(2*pi**2)
	return inft

		

class HardSpheresPY(object):
	def __init__(self, arg):
		self.eta = arg
	# Percus-Yevick Terms  (see for instance D. Henderson Condensed Matter Physics 2009, Vol. 12, No 2, pp. 127-135)
	# Percus-Yevick Terms  (see for instance D. Henderson Condensed Matter Physics 2009, Vol. 12, No 2, pp. 127-135)
	def c0(self):
	    return -(1.+2.*self.eta)**2/(1.-self.eta)**4
	def c1(self):
	    return 6.*self.eta*(1.+self.eta*0.5)**2/(1.-self.eta)**4
	def c3(self):
	    return self.eta*0.5*self.c0()
	 

	def c(self,r):
		"""Hard spheres correlation function c(r) at given packing fraction self.eta.
			SeeM. S. Wertheim Phys. Rev. Lett. 10, 321 (1963)
			"""
		if r>1:
			return 0
		else:
			return self.c0()+self.c1()*r+self.c3()*r**3

	def cc(self,r):
		return vectorize(self.c)(r)
	# Spherical Fourier Transforms (using the liquid isotropicity)

	def PercusYevick(self,dr=0.0005,plot=True,filename="g_of_r.txt", npoints =1024*5):
		# number density
		rho=6./pi*self.eta
		# getting the direct correlation function c(r) from the analytic Percus-Yevick solution
	
		# space discretization dr
		
		r=pl.arange(1,npoints+1,1 )*dr
		# reciprocal space discretization (highest available frequency)
		dk=1/r[-1]
		k=pl.arange(1,npoints+1,1 )*dk
		# direct correlation function c(r)
		c_direct=self.cc(r)
		# getting the Fourier transform
		ft_c_direct=spherical_FT(c_direct, k,r,dr)
		# using the Ornstein-Zernike equation, getting the structure factor
		ft_h=ft_c_direct/(1.-rho*ft_c_direct)
		# inverse Fourier transform
		h=inverse_spherical_FT(ft_h, k,r,dk)
		# print h
		# # radial distribution function
		gg=h+1.0
		# clean the r<1 region
		g=pl.zeros(len(gg))
		g[r>=1]=gg[r>=1]
		# save the cleaned version
		# print (g)
		pl.savetxt(filename, column_stack((r,g)))
		# plots
		if plot:
			pl.plot(r,g)
			pl.ylabel("g(r)")
			pl.xlabel("r")
			# pl.show()
		return r,g

# a = HardSpheresPY(0.5)
# r,g = a.PercusYevick(dr=0.00001,plot=True, npoints = 1024*10)
# pl.show()
