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
N=600
# cubic box side:
L=47.13*mu
# the volume is
V=L**3


import numpy as np

diam=np.random.normal(mu, sigma, size=N)
print len(diam[diam<0])
# random coordinates
x=np.random.uniform(0,L,size=N)
y=np.random.uniform(0,L,size=N)
z=np.random.uniform(0,L,size=N)

packing=np.pi/6*(diam**3).sum()/V

print "    * The packing fraction is", packing

print "    * Writing file",filename,"..."
with open(filename,'w') as fw:
    # write the header
    fw.write("LAMMPS Description\n\n")
    fw.write("%d atoms\n\n"%N)
    # every atom is a type
    fw.write("%d atom types\n"%N)
    fw.write("""
0 %g xlo xhi
0 %g ylo yhi
0 %g zlo zhi
\n"""%(L,L,L))
    # all the atoms have different masses, but equal mass density 1
    fw.write("Masses\n\n")
    for i in xrange(N):
        fw.write("%d %g\n"%(i+1,np.pi/6*diam[i]**3))

    fw.write("\n")
    # specify the interaction coefficient for the Morse potential:
    # in lammps's documentation they are:
    # d0 alpha r0 cutoff
    # in Paddy's notation:
    # epsilon rho0 sigma cutoff
    fw.write("PairIJ Coeffs # morse\n\n")
    for i in xrange(N):
        # fw.write("%d %g %g %g \n"%(i+1,epsilon,rho0,diam[i]))#, r_cut_coeff*diam[i] ))
        for j in xrange(i, N):
            mixed_diam=0.5*(diam[i]+diam[j])
            fw.write("%d %d %g %g %g %g\n"%(i+1,j+1,epsilon,rho0,mixed_diam, r_cut_coeff*mixed_diam ))
    
    fw.write("\nAtoms\n\n")
    # write random coordinates
    for i in xrange(N):
        fw.write("%d %d %g %g %g\n"%(i+1,i+1, x[i],z[i],y[i]))

print "    Compressing..."
import os
# os.system("gzip "+filename)
print "    ...done."

