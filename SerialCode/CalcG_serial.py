from advs_serial import *
from CalcTriStrains_serial import *
import numpy as np
import time


def CalcG(c,v,Obs):
	# some stuff for the transition from MATLAB
	v = v-1 # index from 0
	v = np.array(v,'int') 
	sx = Obs[:,0]
	sy = Obs[:,1]
	if Obs.shape[1]==2: #for tiny dataset observations are on a plane
		sz=np.zeros(len(Obs[:,1]))
	if Obs.shape[1]==3:	
		sz = Obs[:,2]
	
	pr = 0.25
		
	ss = 1
	ts = 0
	ds = 0
	totalkernel=0.0

	totaltransfer=0.0
	G = np.zeros([6*len(sx), 3*len(v)])

	for i in range(len(v)):
	
		x = c[v[i,:],0]
		y = c[v[i,:],1]
		z = c[v[i,:],2]

		[S_xx, S_yy, S_zz, S_xy, S_xz, S_yz, kerneltime] = CalcTriStrains(sx, sy, sz, x, y, z, pr, 1, 0, 0)
		totalkernel=totalkernel+kerneltime
		G[::6,3*i-2]  = S_xx;
		G[1::6,3*i-2] = S_yy;
		G[2::6,3*i-2] = S_zz;
		G[3::6,3*i-2] = S_xy;
		G[4::6,3*i-2] = S_xz;
		G[5::6,3*i-2] = S_yz;

		[S_xx, S_yy, S_zz, S_xy, S_xz, S_yz, kerneltime] = CalcTriStrains(sx, sy, sz, x, y, z, pr, 0, 1, 0)
		totalkernel=totalkernel+kerneltime
		G[::6,3*i-1]  = S_xx;
		G[1::6,3*i-1] = S_yy;
		G[2::6,3*i-1] = S_zz;
		G[3::6,3*i-1] = S_xy;
		G[4::6,3*i-1] = S_xz;
		G[5::6,3*i-1] = S_yz;

		[S_xx, S_yy, S_zz, S_xy, S_xz, S_yz, kerneltime] = CalcTriStrains(sx, sy, sz, x, y, z, pr, 0, 0, 1)
		totalkernel=totalkernel+kerneltime
		G[::6,3*i]  = S_xx;
		G[1::6,3*i] = S_yy;
		G[2::6,3*i] = S_zz;
		G[3::6,3*i] = S_xy;
		G[4::6,3*i] = S_xz;
		G[5::6,3*i] = S_yz;
     
	#filename = "G_parallel_step3_%s" % MATRIX_SIZE
	#np.savetxt(filename, G)
	

	return G, totalkernel
	

	

