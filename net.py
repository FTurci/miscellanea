import numpy as np
import sys
import glob

def read_and_split(filename,N):
    # skip the first line
    with open(filename,'r') as fin:
        fin.readline()
        data=np.genfromtxt(fin,invalid_raise=False,usecols=[0],dtype='S1')
        frames=data.reshape((N+2, len(data)/(N+2)), order='F')
    return frames[2:,:]

def cluster(info):
    return (info=='C')+(info=='D')


if len(sys.argv)<2:
   print "Usage python net.py file.xyz"
   sys.exit(0) 
# input a file name
xyz_name=sys.argv[1]

# Modify this list in order to have your favourite hierarchy
priority_list=['FCC','HCP','10B','7A','6A', '5A']

# read number of particles
with open(xyz_name, 'r') as xyz:
    N=int(xyz.readline())

TCCinfo={}
netTCCinfo={}
gross,net={},{}

for species in priority_list:
    print "Reading", species,"..."
    TCCinfo[species]=read_and_split(glob.glob(xyz_name+"*"+species)[0],N)

# testing 
test=cluster(TCCinfo[priority_list[0]])

netTCCinfo[priority_list[0]]=np.copy(test)
net[priority_list[0]]=netTCCinfo[priority_list[0]].sum(axis=0)/float(N)
gross[priority_list[0]]=net[priority_list[0]]


for species in priority_list[1:]:
    c=cluster(TCCinfo[species])
    netTCCinfo[species]=np.logical_not(test)*c
    test+=c
    net[species]=netTCCinfo[species].sum(axis=0)/float(N)
    gross[species]=c.sum(axis=0)/float(N)

for species in priority_list:
    np.save("net_"+species,netTCCinfo[species])

# save the rest
np.save("net_other", np.logical_not(test))

# read xyz data
nframes=netTCCinfo[priority_list[0]].shape[1]
xyz=np.genfromtxt(xyz_name,skiprows=2, invalid_raise=False,usecols=[1,2,3])
xyz=xyz.reshape([nframes, N,3])
# write a single xyz file
with open("merged.txt",'w') as fout:
    for frame in range():
        fout.write("%d\nAtoms\n")
        for p in range(N):
            found=False
            iterator=iter(priority_list)
            while found==False:
                species=iterator.next()
                if netTCCinfo[species][p,frame]:
                    fout.write(species+" %g %g %g\n"%(xyz[frame, p,0],xyz[frame, p,1],xyz[frame, p,2]))


with open(xyz_name+"_net",'w') as fout:
    fout.write("Species\tGross\tNet\n")
    for species in priority_list:
        fout.write(species+":\t%g\t%g\n"%(gross[species], net[species]))

