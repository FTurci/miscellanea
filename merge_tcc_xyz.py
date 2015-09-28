import sys 
import numpy as np

if (len(sys.argv)<3):
    print "Usage: python merge_columns.py your_coords.xyz your_raw_file_raw"
    sys.exit(0)

XYZ=np.genfromtxt(sys.argv[1], invalid_raise=False, skiprows=2, usecols=[1,2,3])
with open(sys.argv[1], 'r') as fXYZ:
    N=int(fXYZ.readline())

with open(sys.argv[2], 'r') as fRAW:
    # skip the first line
    fRAW.readline()
    # read the rest
    RAW_=fRAW.readlines()
    # trim only the types information: skipping two lines per frame
    nFrames=len(RAW_)/(N+2)
    print nFrames
    RAW=[]
    for frame in range(nFrames):
        for p in range(N):
            RAW.append(RAW_[frame*(N+2)+p+2][0])

print "I am writing",nFrames," frames"


with open("merged_"+sys.argv[1], 'w') as fw:
    for frame in range(nFrames):
        fw.write("%d\nAtoms\n"%N)
        for p in range(N):
            i=frame*N+p
            fw.write("%s %g %g %g\n"%(RAW[i],XYZ[i,0],XYZ[i,1],XYZ[i,2]))

