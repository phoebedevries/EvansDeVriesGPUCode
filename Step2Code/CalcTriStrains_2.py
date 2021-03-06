from advs_parallel_2 import *
import numpy as np
import time
def CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds):

	##################################################
	#
	### Inputs
	#
	# sx : x-coordinates of observation points
	# sy : y-coordinates of observation points
		# sz : z-coordinates of observation points
	# x  : x-coordinates of triangle vertices 
	# y  : y-coordinates of triangle vertices
	# z  : z-coordinates of triangle vertices
	# pr : Poisson's ratio
	# ss  : strike slip displacement
	# ts  : tensile slip displacement
	# ds  : dip slip dislplacement
	#
	#################################################

	nsta = len(sx)

	# Calculate the slip vector in XYZ coordinates

	normVec = np.cross(np.array([x[1],y[1],z[1]])-np.array([x[0],y[0],z[0]]),np.array([x[2],y[2],z[2]])-np.array([x[0],y[0],z[0]]))
	normVec = normVec/np.linalg.norm(normVec)
        kerneltime=0.0
        transfertime=0.0
	# Enforce clockwise circulation
	if normVec[2] < 0:
		normVec = -normVec
		x[1], x[2] = x[2], x[1]
		y[1], y[2] = y[2], y[1]
		z[1], z[2] = z[2], z[1]

	strikeVec = np.array([-np.sin(np.arctan2(normVec[1],normVec[0])), np.cos(np.arctan2(normVec[1],normVec[0])), 0.0])
	dipVec = np.cross(normVec,strikeVec)
	slipComp = np.array([ss, ds, ts])
	slipVec = np.dot(np.array([strikeVec, dipVec, normVec]).T,slipComp)

	# Solution Vectors
	S_xx = np.zeros(len(sx))
	S_yy = np.zeros(len(sx))
	S_zz = np.zeros(len(sx))
	S_xy = np.zeros(len(sx))
	S_xz = np.zeros(len(sx))
	S_yz = np.zeros(len(sx))

	# Add a copy of the first vertex to the end of the vertex list for indexing
	x = np.append(x,x[0])
	y = np.append(y,y[0])
	z = np.append(z,z[0])

	pi = np.pi

	# Calculate dislocations
	for iTri in range(3):
		# Calculate strike and dip of current leg
		strike = (180/pi)*( np.arctan2( y[iTri+1]-y[iTri] , x[iTri+1]-x[iTri] ) )
		segMapLength = np.sqrt( (x[iTri]-x[iTri+1])**2 + (y[iTri]-y[iTri+1])**2 )
		# Rotate by strike
		alpha = (pi/180)*(-strike);
		rx = np.cos(alpha)*( x[iTri+1] - x[iTri] ) - np.sin(alpha)*( y[iTri+1] - y[iTri] )
		ry = np.sin(alpha)*( x[iTri+1] - x[iTri] ) + np.cos(alpha)*( y[iTri+1] - y[iTri] )
	
		dip = (180/pi)*(np.arctan2(z[iTri+1] - z[iTri], rx))

		if dip >= 0: 
			beta = (pi/180)*(90-dip);
			if beta > pi/2:
				beta = (pi/2)-abs(beta)
		else:
			beta = (-pi/180)*(90+dip)
			if beta < -pi/2:
				beta = (pi/2) - abs(beta)
	
		ssVec = np.array([np.cos(strike/180*pi), np.sin(strike/180*pi), 0])
		tsVec = np.array([-np.sin(strike/180*pi), np.cos(strike/180*pi), 0])
		dsVec = np.cross(ssVec, tsVec)
		lss = np.dot(slipVec,ssVec)
		lts = np.dot(slipVec,tsVec)
		lds = np.dot(slipVec,dsVec)

		if (abs(beta) > 0.000001) and (abs(beta - pi) > 0.000001):
			# First angular dislocation
			sx1 = np.cos(alpha)*( sx - x[iTri] ) - np.sin(alpha)*( sy - y[iTri] )
			sy1 = np.sin(alpha)*( sx - x[iTri] ) + np.cos(alpha)*( sy - y[iTri] )
			sz1 = sz-z[iTri]
			a1 = np.array(z[iTri])

			# Second angular dislocation
			sx2 = np.cos(alpha)*( sx - x[iTri+1] ) - np.sin(alpha)*( sy - y[iTri+1] )
			sy2 = np.sin(alpha)*( sx - x[iTri+1] ) + np.cos(alpha)*( sy - y[iTri+1] )
			sz2 = sz-z[iTri+1]
			a2 = np.array(z[iTri+1])

			sx3 = np.append(sx1,sx2);
			sy3 = np.append(sy1,sy2);
			sz3 = np.append(sz1,sz2);
			a3 = np.append(np.tile(a1,len(sx)),np.tile(a2,len(sx)))
			beta3 = np.tile(beta,2*len(sx))
			lss3 = np.tile(lss,2*len(sx))			
			lts3 = np.tile(lts,2*len(sx))			
			lds3 = np.tile(lds,2*len(sx))			
			
			starttime=time.time()
			[c11, c22, c33, c12, c13, c23, ktime] = advs(sx3, sy3, sz3, a3, beta3, pr, lss3, lts3, lds3)
			endtime=time.time()
                        transfertime=transfertime+(endtime-starttime)
                        kerneltime=kerneltime+ktime;
			a11 = c11[0:nsta]
			b11 = c11[nsta::]
			a22 = c22[0:nsta]
			b22 = c22[nsta::]
			a33 = c33[0:nsta]
			b33 = c33[nsta::]
			a12 = c12[0:nsta]
			b12 = c12[nsta::]
			a13 = c13[0:nsta]
			b13 = c13[nsta::]
			a23 = c23[0:nsta]
			b23 = c23[nsta::]

			# Rotate tensors to correct for strike
			bxx = a11-b11
			byy = a22-b22
			bzz = a33-b33
			bxy = a12-b12
			bxz = a13-b13
			byz = a23-b23

			g = (pi/180)*strike
			e11n = (np.cos(g)*bxx - np.sin(g)*bxy)*np.cos(g) - (np.cos(g)*bxy - np.sin(g)*byy)*np.sin(g)
			e12n = (np.cos(g)*bxx - np.sin(g)*bxy)*np.sin(g) + (np.cos(g)*bxy - np.sin(g)*byy)*np.cos(g)
			e13n = np.cos(g)*bxz - np.sin(g)*byz
			e22n = (np.sin(g)*bxx + np.cos(g)*bxy)*np.sin(g) + (np.sin(g)*bxy + np.cos(g)*byy)*np.cos(g)
			e23n = np.sin(g)*bxz + np.cos(g)*byz
			e33n = bzz
		
			# Add the strains from the current leg
			S_xx = S_xx + e11n
			S_yy = S_yy + e22n
			S_zz = S_zz + e33n
			S_xy = S_xy + e12n
			S_xz = S_xz + e13n
			S_yz = S_yz + e23n

	return (S_xx, S_yy, S_zz, S_xy, S_xz, S_yz, kerneltime, transfertime)	
