#/usr/bin/env python3
import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
from numba import njit,prange
plt.style.use('ggplot')

# parameters
k_B = 1
trajectory_file = "traj_unwrapped.xyz"



# notation
# l,l' -> refers to the unit cell in the supercell
# k,k' -> runs over all Natoms
# \alpha, \beta -> x,y,z

cell = np.ones((3,3))
# number of repeats
l = [5,5,5]
# origin of l-th unit cell relative to l=[0,0,0]
r_l = np.dot(l,cell)


def read_xyz(pos=trajectory_file):
	""" Reads .xyz file and create a frame-indexed trajectory array."""
	with open(pos,"r") as f:
		Natoms,Nframes = 0,0
		lines = f.readlines()
		frame = 0
		k = 0
		counter = 0
		for idx,ln in enumerate(lines):
			if counter==0:
				Natoms = int(ln)
				Nframes = int(len(lines)/(Natoms+2))
				traj = np.zeros((Nframes,Natoms,4))
				counter += 1
			elif counter==1:
				k = 0 # atom idx
				counter +=1 
			elif counter>1:
				b = [float(l) for l in ln.split()]
				traj[frame,k,:] = [float(l) for l in ln.split()]
				k += 1
				counter += 1
				if counter==(Natoms+2):
					counter = 0
					frame += 1
	return traj,Natoms,Nframes


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


# @njit()
def greens_func(traj):
	"""	Takes the Fourier transform of the absolute positions for a given vector.
	Averages all frames and calculates the FT Green's function coeffs at each 
	wave vector q."""
	# Nframes = G.shape[0]
	# Natoms = G.shape[1]
	# R_ka = np.mean(traj,axis=0) # average over all frames

	# Natoms = R_ka.shape[0]
	G_ft = np.zeros((Natoms,3,Natoms,3),dtype='complex128') # ka, k'b
	# print("G shape",G_ft.shape)

	# Rsum 
	Rsum = np.mean(traj,axis=0)

	# Rqsum
	Rq = np.fft.fftn(traj)

	# print("Rq shape",Rq.shape)

	# add 2nd term to G
	for k1 in prange(Natoms):
		for k2 in prange(Natoms):
			for a in prange(3):
				for b in prange(3):
					G_ft[k1,a,k2,b] -= np.mean(Rq[:,k1,a]*np.conj(Rq[:,k2,b]),axis=0)


	# sum(R.R)(q)
	# sumRkRk = np.zeros(G_ft.shape)
	# for k1 in range(Natoms):
	# 	for k2 in range(Natoms):
	# 		for a in range(3):
	# 			for b in range(3):
	# 				sumRkRk[k1,a,k2,b] = np.dot(Rsum[k1,a],np.conj(Rsum[k2,b]))
	# sumRqRq = np.fft.fftn(sumRkRk)
	# add first term to G
	# G_ft += sumRqRq


	# shouldn't be like this?
	Rqsum = np.fft.fftn(Rsum)
	# print np.shape(Rsum), np.shape(Rq), np.shape(Rqsum), np.shape(traj), np.shape(Rqsum)
	# sumRkRk = np.zeros(G_ft.shape)
	for k1 in prange(Natoms):
		for k2 in prange(Natoms):
			for a in prange(3):
				for b in prange(3):
					G_ft[k1,a,k2,b] +=  Rqsum[k1,a]*np.conj(Rqsum[k2,b])
	# sumRqRq = np.fft.fftn(sumRkRk)





	# to be removed
	# R_ka_kb = R_ka*R_kb
	# R_ka_ft = np.fft.fft(R_ka)
	# R_kb_ft = np.conj(np.fft.fft(R_kb))
	# R_ka_kb_ft = np.fft.fft(R_ka*R_kb)
	# G_ft = np.zeros((Natoms,3,Natoms,3)) # ka, k'b
	# G_ft = R_ka_kb_ft-R_ka_ft*R_kb_ft
	print(G_ft.shape)
	return G_ft

# @njit()
def force_constants(G,T=20):
	""" Calculates force constants $\Phi_{lk\alpha,l'k'\beta}$ """
	# phi = np.zeros(np.shape(G))

	# check if G is hermitian 
	# !! PROBLEM should check for all atoms separately
	# if (np.conj(G)==G).all(): # check if G is hermitian 
		# print("Matrix is Hermitian")
	# else:
		# print("Matrix is NOT Hermitian\n",np.conj(G)==G)
	# 	phi = k_B * T* G
	# else:

	# for k1 in range(Natoms):
	# 	for k2 in range(Natoms):
	# 		for a in range(3):
	# 			for b in range(3):
	# 				print("G shape atom_i ={} atom_j={} |".format(k1,k2),G[k1,:,k2,:].shape)
	# 				print("== ==\n",G[k1,:,k2,:],np.conj(G[k1,:,k2,:].T))
	# 		phi[k1,:,k2,:] = k_B * T* np.linalg.inv(G[k1,:,k2,:])
	return G


def eigenfreqs(traj,M=1.):
	G_ft = greens_func(traj)
	# print(G_ft.shape)
	# phi_ft = force_constants(G_ft)
	phi_ft = force_constants(G_ft)
	D = np.zeros(phi_ft.shape)
	D = 1/np.sqrt(M)* phi_ft
	omega_sq = np.zeros((Natoms,3),dtype='float64')
	for k1 in prange(Natoms):
		for k2 in prange(Natoms):
			# print(D[k1,:,k2,:].shape)
			# print(np.linalg.eig(D[k1,:,k2,:]))
			eigenvals,eigenvecs = np.linalg.eigh(D[k1,:,k2,:])
			print("== EIGENVALUES ==\n",eigenvals)
			omega_sq[k1,:] = eigenvals
	return np.sqrt(omega_sq)


def main():
	global Natoms
	global Nframes
	a = 1. # lattice constant
	# high symmetry points for fcc Gamma -> X -> W -> K -> Gamma -> L
	gamma = np.array([0,0,0])
	X = np.array([0,2*np.pi/a,0])
	W = np.array([np.pi/a,2*np.pi/a,0])
	K =  np.array([3*np.pi/(2*a),3*np.pi/(2*a),0])
	L = np.array([np.pi/a,np.pi/a,np.pi/a])	
	pt = highsymm_path(np.array([gamma,X,W,K,gamma,L]))
	# plt.plot(pt,"-.")
	# plt.legend(["1/x","1/y","1/z"])
	# plt.show()



	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	traj,Natoms,Nframes = read_xyz()
	traj = traj[:,:,1:]

	# Gf=greens_func(traj)
	freqs = eigenfreqs(traj)
	print(" == FREQUENCIES (omega(q)) ==\n",freqs)
	# ax.scatter(freqs[:,0],freqs[:,1],freqs[:,2])
	# plt.show()

	return None

if __name__ == '__main__':
	main()