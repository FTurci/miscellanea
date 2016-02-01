# October 2015
# 
# Writing a LAMMPS input file for a system with Morse interaction and gaussian polydispersity

filename="polydisperse_morse_input.lmp"
epsilon=1
rho0=33.0
r_cut_coeff=1.4
mu=1.0
sigma=0.04
polydispersity=sigma/mu
print "    * The imposed polydispersity is", polydispersity
# number of atoms
N=10500
# cubic box side:
L=47.13*mu
# the volume is
V=L**3


import numpy as np
from scipy.stats import norm


Ds=np.arange(mu-3*sigma, mu+3*sigma+sigma, sigma)-sigma/2.
diams=Ds[:-1]+(Ds[1]-Ds[0])/2
label_types=np.arange(len(diams)+1)+1
cumulative=norm.cdf(Ds,mu,sigma)
cumulative[0]=0
cumulative[-1]=1
Types=[]
for p in range(N):
    random_num=np.random.uniform(0,1)
    Test=cumulative<random_num
    Types.append(label_types[Test][-1])


Num_Types=len(Ds)-1

# random coordinates
x=np.random.uniform(0,L,size=N)
y=np.random.uniform(0,L,size=N)
z=np.random.uniform(0,L,size=N)

packing=np.pi/6*(diams**3).sum()/V

print "    * The packing fraction is", packing

print "    * Writing file",filename,"..."
with open(filename,'w') as fw:
    # write the header
    fw.write("LAMMPS Description\n\n")
    fw.write("%d atoms\n\n"%N)
    # every atom is a type
    fw.write("%d atom types\n"%Num_Types)
    fw.write("""
0 %g xlo xhi
0 %g ylo yhi
0 %g zlo zhi
\n"""%(L,L,L))
    # all the atoms have different masses, but equal mass density 1
    fw.write("Masses\n\n")
    for i in xrange(Num_Types):
        fw.write("%d %g\n"%(i+1,np.pi/6*diams[i]**3))

    fw.write("\n")
    # specify the interaction coefficient for the Morse potential:
    # in lammps's documentation they are:
    # d0 alpha r0 cutoff
    # in Paddy's notation:
    # epsilon rho0 sigma cutoff
    fw.write("PairIJ Coeffs # morse\n\n")
    for i in xrange(Num_Types):
        # fw.write("%d %g %g %g \n"%(i+1,epsilon,rho0,diam[i]))#, r_cut_coeff*diam[i] ))
        for j in xrange(i, Num_Types):
            diam_i=diams[i]
            diam_j=diams[j]
            mixed_diam=0.5*(diam_i+diam_j)
            fw.write("%d %d %g %g %g %g\n"%(i+1,j+1,epsilon,rho0,mixed_diam, r_cut_coeff*mixed_diam ))
    
    fw.write("\nAtoms\n\n")
    # write random coordinates
    for i in xrange(N):
        fw.write("%d %d %g %g %g\n"%(i+1,Types[i], x[i],z[i],y[i]))

print "    Compressing..."
import os
# os.system("gzip "+filename)
print "    ...done."

