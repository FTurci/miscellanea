import pylab as pl
from numpy import pi

# Percus-Yevick Terms  (see for instance D. Henderson Condensed Matter Physics 2009, Vol. 12, No 2, pp. 127-135)
def c0(eta):
    return -(1.+2.*eta)**2/(1.-eta)**4
def c1(eta):
    return 6.*eta*(1.+eta*0.5)**2/(1.-eta)**4
def c3(eta):
    return eta*0.5*c0(eta)
def cc(r,eta):
    if r>1:
        return 0
    else:
        return c0(eta)+c1(eta)*r+c3(eta)*r**3
# Spherical Fourier Transforms (using the liquid isotropicity)
def spherical_FT(f,k,r,dr):
    ft=pl.zeros(len(k))
    for i in range(len(k)):
        ft[i]=4.*pi*pl.sum(r*pl.sin(k[i]*r)*f*dr)/k[i]
    return ft

def inverse_spherical_FT(ff,k,r,dk):
    ift=pl.zeros(len(r))
    for i in range(len(r)):
        ift[i]=pl.sum(k*pl.sin(k*r[i])*ff*dk)/r[i]/(2*pi**2)
    return ift

def PercusYevickHS(phi,plot=True,filename="g_of_r.txt", factor=8, shift=0):
    # number density
    rho=6./pi*phi
    # getting the direct correlation function c(r) from the analytic Percus-Yevick solution
    # vectorizing the function
    c=pl.vectorize(cc)
    # space discretization
    dr=0.0005
    r=pl.arange(1,1024*factor+1,1 )*dr
    # reciprocal space discretization (highest available frequency)
    dk=1/r[-1]
    k=pl.arange(1,1024*factor+1,1 )*dk
    # direct correlation function c(r)
    c_direct=c(r,phi)
    # getting the Fourier transform
    ft_c_direct=spherical_FT(c_direct, k,r,dr)
    # using the Ornstein-Zernike equation, getting the structure factor
    ft_h=ft_c_direct/(1.-rho*ft_c_direct)
    # inverse Fourier transform
    h=inverse_spherical_FT(ft_h, k,r,dk)
    # print h
    # # radial distribution function
    gg=h+1
    # clean the r<1 region
    g=pl.zeros(len(gg))
    g[r>=1]=gg[r>=1]
    # save the cleaned version
    pl.savetxt(filename, list(zip(r,g)))
    # plots
    if plot:
        pl.plot(r,g+shift)
        pl.ylabel("g(r)")
        pl.xlabel(r"$r/\sigma$")

# call the function
# PercusYevickHS(0.2)
# PercusYevickHS(0.3)
# PercusYevickHS(0.45)
# PercusYevickHS(0.525)
# PercusYevickHS(0.53)
# PercusYevickHS(0.54)
# PercusYevickHS(0.55)
for i,phi in enumerate(pl.linspace(0.1,0.55,3)):
    PercusYevickHS(phi, factor=12,shift=i)
# pl.xlim(1,2)
pl.savefig("profs.svg")
pl.show()
