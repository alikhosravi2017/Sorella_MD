#/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
plt.style.use('ggplot')

# parameters
k_B = 1


# notation
# l,l' -> refers to the unit cell in the supercell
# k,k' -> runs over all Natoms
# \alpha, \beta -> x,y,z

cell = np.ones((3,3))
# number of repeats
l = [5,5,5]
# origin of l-th unit cell relative to l=[0,0,0]
r_l = np.dot(l,cell)

@njit()
def force_constants(G,T):
	""" Calculates force constants $\Phi_{lk\alpha,l'k'\beta}$ """
	phi = k_B * T* np.linalg.inv(G)
	return phi

# @njit()
def equidist(p1, p2, npoints=20):
	""" Creates an array of equidistant points between two N-dim points."""
	temp = np.zeros((npoints,p1.shape[0]))
	for i in range(p1.shape[0]):
		temp[:,i] = np.linspace(p1[i],p2[i],npoints)
	return temp

# @njit()
def highsymm_path(symm_points,npoints=200):
	""" Generates high symmetry path 
	along the given points.
	a -> lattice constant"""
	# symm_points = np.array(symm_points)
	path = np.zeros((npoints*(symm_points.shape[0]-1),symm_points.shape[1]))

	for i in range(len(symm_points)-1):
		temp = equidist(symm_points[i],symm_points[i+1],npoints)
		# print(temp.shape)
		z = 0
		for j in range(i*npoints,(i+1)*npoints):
			path[j,:] = temp[z,:]
			z += 1

	# hardcoded (to be removed)
	# gammaX = equidist(gamma,X,20)
	# XW = equidist(X,W,20)
	# WK = equidist(W,K,20)
	# Kgamma = equidist(K,gamma,20)
	# gammaL = equidist(gamma,L,20+1)
	# path = np.vstack((gammaX,XW,WK,Kgamma,gammaL))

	return path

a = 1. # lattice constant
# high symmetry points for fcc Gamma -> X -> W -> K -> Gamma -> L
gamma = np.array([0,0,0])
X = np.array([0,2*np.pi/a,0])
W = np.array([np.pi/a,2*np.pi/a,0])
K =  np.array([3*np.pi/(2*a),3*np.pi/(2*a),0])
L = np.array([np.pi/a,np.pi/a,np.pi/a])	
pt = highsymm_path(np.array([gamma,X,W,K,gamma,L]))
plt.plot(pt,"-.")
plt.legend(["1/x","1/y","1/z"])
plt.show()
