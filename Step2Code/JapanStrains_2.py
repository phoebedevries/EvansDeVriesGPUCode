from advs_parallel_2 import *
from CalcTriStrains_2 import *
from CalcG_2 import *
import numpy as np
import time


PROBLEM_SIZE = "TINY"

totalstart = time.time()
if PROBLEM_SIZE=="TINY":
	# Set up inputs for script instead of function (function later)
	c = np.loadtxt('c_tiny.txt') # coordinates of triangle vertices
	v = np.loadtxt('v_tiny.txt') # indices of triangle vertices from MATLAB (indexing from 1)
	slip = np.loadtxt('Slip_small.txt')
	Obs = np.loadtxt('Obs_tiny.txt') # Observation locations

elif PROBLEM_SIZE=="SMALL":
	# Set up inputs for script instead of function (function later)
	c = np.loadtxt('c_small.txt') # coordinates of triangle vertices
	v = np.loadtxt('v_small.txt') # indices of triangle vertices from MATLAB (indexing from 1)
	slip = np.loadtxt('Slip_small.txt')
	Obs = np.loadtxt('Obs_small.txt') # Observation locations

elif PROBLEM_SIZE=="MEDIUM":
	# Set up inputs for script instead of function (function later)
	c = np.loadtxt('c_med.txt') # coordinates of triangle vertices
	v = np.loadtxt('v_med.txt') # indices of triangle vertices from MATLAB (indexing from 1)
	slip = np.loadtxt('Slip_med.txt')
	Obs = np.loadtxt('Obs_med.txt') # Observation locations

elif PROBLEM_SIZE=="FINAL":
	# Set up inputs for script instead of function (function later)
	c = np.loadtxt('c_final.txt') # coordinates of triangle vertices
	v = np.loadtxt('v_final.txt') # indices of triangle vertices from MATLAB (indexing from 1)
	slip = np.loadtxt('Slip_final.txt')
	Obs = np.loadtxt('Obs_med.txt') # Observation locations

# Calculate G
[G,totalkernel,totaltransfer] = CalcG(c,v,Obs)

# Calculate Strains
slipfull = np.zeros(3*len(slip))
slipfull[1::3] = slip

E = np.dot(G,slipfull)
Exx = E[0::6]
Eyy = E[1::6]
Ezz = E[2::6]
Exy = E[3::6]
Exz = E[4::6]
Eyz = E[5::6]

exxfilename = "Exx_parallel_step2_%s" % PROBLEM_SIZE
eyyfilename = "Eyy_parallel_step2_%s" % PROBLEM_SIZE
ezzfilename = "Ezz_parallel_step2_%s" % PROBLEM_SIZE
exyfilename = "Exy_parallel_step2_%s" % PROBLEM_SIZE
exzfilename = "Exz_parallel_step2_%s" % PROBLEM_SIZE
eyzfilename = "Eyx_parallel_step2_%s" % PROBLEM_SIZE

totalend=time.time()   

print "\ntotal time to execute advs kernel: %f" % totalkernel
print "total time to execute kernel including setup and transfer : %f" % totaltransfer
print "total time execute the entire program and build G : %f\n" % (totalend-totalstart)
print "Please see the files for the outputed strains E:\n %s\n %s\n %s\n %s\n %s\n %s\n" % (exxfilename, eyyfilename, ezzfilename, exyfilename, exzfilename, eyzfilename)

np.savetxt(exxfilename, Exx)
np.savetxt(eyyfilename, Eyy)
np.savetxt(ezzfilename, Ezz)
np.savetxt(exyfilename, Exy)
np.savetxt(exzfilename, Exz)
np.savetxt(eyzfilename, Eyz)
	

# Yay!
