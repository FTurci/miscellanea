import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad 
def lennard_jones(r,eps=1,sigma=1):
    return 4*eps*((sigma/r)**12-(sigma/r)**6)

def wca(r):
    rcut = 2**(1/6)
    if r>rcut:
        return 0
    else:
        return lennard_jones(r)+1

v_wca = np.vectorize(wca)


def baker_henderson(T):

    beta = 1./T
    def integrand(r):
        return 1-np.exp(-beta*wca(r))

    return quad(integrand, 0,50)
# plt.plot(r,wca(r))
# plt.show()

print("Baker-Henderson diameter",baker_henderson(1.0))