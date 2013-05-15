import pylab, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import sys
#membrane-0.50-2.50-2.50-3
#f_membrane-%4.02f-%4.02f-%4.02f-%d.dat", A,B,C,rand_seed
#print sys.argv[1]

fname = "shape-" + str(sys.argv[1]) + "/f_membrane-" + str(sys.argv[2]) + "-" + str(sys.argv[3]) + "-" + str(sys.argv[4]) + "-" + str(sys.argv[5]) + ".dat"
print fname
A = np.loadtxt(fname, dtype = float)
width = int(sqrt(len(A)))
print width
Y = reshape(A[:,0],(width,width))
Z = reshape(A[:,1],(width,width))
f = reshape(A[:,2],(width,width))

contour(Y,Z,f)
colorbar()
show()
