#/usr/bin/env python3
import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
from numba import njit,prange
plt.style.use('ggplot')

# parameters
k_B = 1
T =1
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
		traj_built_flag = True
		for idx,ln in enumerate(lines):
			if counter==0:
				if traj_built_flag :
					Natoms = int(ln)
					Nframes = int(len(lines)/(Natoms+2))
					traj = np.zeros((Nframes,Natoms,4))
					traj_built_flag = False
				counter += 1
			elif counter==1:
				k = 0 # atom idx
				counter +=1 
			elif counter>1:
				b = [float(l) for l in ln.split()]
				# print [float(l) for l in ln.split()]
				traj[frame,k,:] = [float(l) for l in ln.split()]
				# print traj[frame,k,:]
				k += 1
				counter += 1
				if counter==(Natoms+2):
					# print Natoms, ln
					counter = 0
					frame += 1
	return traj,Natoms,Nframes


# @njit()
def equidist(p1, p2, npoints=20):
	""" Creates an array of equidistant points between two N-dim points."""
	temp = np.zeros((npoints,p1.shape[0]))
	# loop over x,y,z dimensions
	for i in range(p1.shape[0]):
		temp[:,i] = np.linspace(p1[i],p2[i],npoints)
	return temp

# @njit()
def highsymm_path(symm_points,npoints=100):
	""" Generates high symmetry path 
	along the given points."""

	path = np.zeros((npoints*(symm_points.shape[0]-1),symm_points.shape[1]))

	# loop over each high symm point
	for i in range(len(symm_points)-1):
		temp = equidist(symm_points[i],symm_points[i+1],npoints)
		# print(temp.shape)
		z = 0
		# loop over each sample point along the path
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


def plot_disp(bands):
	""" Plots phonon dispersion.
	bands shape: No. bands x No. path points x No. dimensions"""
	x = np.arange(0,bands.shape[1])
	for band_i in range(bands.shape[0]):
		en = bands[band_i,:,0]
		plt.plot(x,en)
	plt.show()
#	plt.savefig('phonon_dispersion.pdf')
	return None


def FT(POS,pt):
	'''
	Calculate Fourier transform of position. ONLY for 1 FRAME
	:param traj: in this dimension: POS[natom,3]
	:param pt: high symmetry path
	:return:
	'''
	# print "POS",POS
	kir = np.zeros((len(pt),3),dtype = "complex_")
	if len(POS.shape) == 2:
		for ii in range(len(pt)):
			qq = pt[ii]
			exponential = np.exp(-1j * np.sum(qq * POS, axis=1))
			# print  np.sum(POS[:,0] * exponential)
			kir[ii,0] = Natoms_root_rev*np.sum( POS[:,0] * exponential )
			kir[ii,1] = Natoms_root_rev*np.sum( POS[:,1] * exponential )
			kir[ii,2] = Natoms_root_rev*np.sum( POS[:,2] * exponential )
	else:   raise ValueError
	return kir


# @njit()
def greens_func(traj,pt):
	"""	Takes the Fourier transform of the absolute positions for a given vector.
	Averages all frames and calculates the FT Green's function coeffs at each 
	wave vector q."""
	# Nframes = G.shape[0]
	# Natoms = G.shape[1]
	# R_ka = np.mean(traj,axis=0) # average over all frames

	G_ft = np.zeros((len(pt),3,3),dtype='complex128') # ka, k'b
	# print("G shape",G_ft.shape)
	# print 'traj is',traj
	# For first term
	for fram in range(Nframes):
		Rq      =  FT(traj[fram],pt)
		Rq_star =  np.conj(Rq)
		for qq in range(len(pt)):
			for alpha in range(3):
				for beta in range(3):
					G_ft[qq,alpha,beta]+= Rq[qq,alpha]*Rq_star[qq,beta]
					# print G_ft[qq]
	G_ft*(1.0/Nframes)

	# For Second term
	R_mean = np.mean(traj,axis=0)
	R_mean_q = FT(R_mean,pt)
	R_mean_q_star = np.conj(R_mean_q)
	for qq in range(len(pt)):
		for alpha in range(3):
			for beta in range(3):
				G_ft[qq,alpha, beta] -= R_mean_q[qq, alpha] * R_mean_q_star[qq, beta]


	print("G_ft.shape=",G_ft.shape)
	return G_ft

# @njit()
def force_constants(G):
	""" Calculates force constants $\Phi_{lk\alpha,l'k'\beta}$ """
	# phi = np.zeros(np.shape(G))

	# check if G is hermitian 
	# !! PROBLEM should check for all atoms separately
	for qq in range(G.shape[0]):
		if (np.round(np.conj(G[qq]),1)==np.round(G[qq],1)).all(): # check if G is hermitian
			print(G[qq])
			print("Matrix is Hermitian, and Determinant is=",np.linalg.det(G[qq]))
		else:
			# print("Matrix is NOT Hermitian\n",np.conj(G)==G)
			print("Matrix is NOT Hermitian")
			# print "G.conj is :\n",np.conj(G[qq])
			print("G is :\n",G[qq][0])
			print("G is :\n",G[qq][1])
			print("G is :\n",G[qq][2])
			#exit()
	# 	phi = k_B * T* G
	# else:

	# for k1 in range(Natoms):
	# 	for k2 in range(Natoms):
	# 		for a in range(3):
	# 			for b in range(3):
	# 				print("G shape atom_i ={} atom_j={} |".format(k1,k2),G[k1,:,k2,:].shape)
	# 				print("== ==\n",G[k1,:,k2,:],np.conj(G[k1,:,k2,:].T))
	# 		phi[k1,:,k2,:] = k_B * T* np.linalg.inv(G[k1,:,k2,:])

	#### FROM EQ. 17 of the phonons paper
	Phi = np.zeros(G.shape) # ka, k'b
	for qq in range(G.shape[0]):
		Phi[qq] = k_B*T *np.linalg.inv(G[qq])
	####
	return Phi


def eigenfreqs(traj,pt,M=1.):
	G_ft = greens_func(traj,pt)
	# print(G_ft.shape)
	phi_ft = force_constants(G_ft)
	# D = np.zeros(phi_ft.shape)
	D = 1/np.sqrt(M*M)* phi_ft

	omega_sq = np.zeros((len(pt),3),dtype='float64')
	for qq in prange(len(pt)):
		eigenvals,eigenvecs = np.linalg.eigh(D[qq])
		print("== EIGENVALUES ==\n",eigenvals)
		omega_sq[qq] = eigenvals
	print("succeed")
	return np.sqrt(omega_sq)


def main():
	global Natoms
	global Natoms_root_rev
	global Nframes
	a = 1. # lattice constant
	# high symmetry points for fcc Gamma -> X -> W -> K -> Gamma -> L
	gamma = np.array([0,0,0])
	X = np.array([0,2*np.pi/a,0])
	W = np.array([np.pi/a,2*np.pi/a,0])
	K =  np.array([3*np.pi/(2*a),3*np.pi/(2*a),0])
	L = np.array([np.pi/a,np.pi/a,np.pi/a])	
	U = np.array([np.pi/(2*a),2*np.pi/a,np.pi/(2*a)])
	pt = highsymm_path(np.array([K,gamma,L,W,X,U,X,gamma]))
	# plt.plot(pt,"-.")
	# plt.legend(["1/x","1/y","1/z"])
	# plt.show()


	traj,Natoms,Nframes = read_xyz()
	Natoms_root_rev = 1.0/np.sqrt(Natoms)

	# print "this is traj:",traj[:,:,1:]
	# exit()
	traj = traj[:,:,1:]


	# EXPECTED WORKFLOW
	# inputs: 
	#	trajectory (variable: traj, shape: No. frames x No. atoms x 3)
	# 	high symmetry path (variable: pt, shape: No. path points x 3)
	# output: 
	#	eigenfrequencies (shape: No. eigenfrequencies x No. path points x 3)

	# plot_disp(eigenfrequencies

	# this needs to be changed 
	freqs = eigenfreqs(traj,pt) 
	print(" == FREQUENCIES (omega(q)) ==\n",freqs)
	kline=np.arange(freqs.shape[0])
	plt.plot(kline,freqs[:,0],kline,freqs[:,1],kline,freqs[:,2])
	plt.show()
	#plt.savefig("test.pdf")


	return None

if __name__ == '__main__':
	main()