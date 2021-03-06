#!/usr/bin/env python3
import numpy as np
import time
import matplotlib.pyplot as plt
from numba import njit,prange
plt.style.use('ggplot')

# parameters
folder_path='./'
save_flag= ''
trajectory_file = "traj_unwrapped.xyz"
# folder_path='with_lammps/'
# save_flag= '_10K_1fs_4_4_4'
# trajectory_file = "traj_lammps_10K_1fs.xyz"


sigma_true = 3.4e-10  # m ## This is necessarily to get rid of lj units
kB_true = 1.38064852e-23  #m2 kg s-2 K-1
T = 10
mass = 6.6335209e-26  #kg



@njit(parallel=True,fastmath=True)
def mean(arr):
	summation = np.zeros((arr.shape[1],arr.shape[2]),dtype=np.complex128)
	for j in prange(arr.shape[1]):
		for frame in range(arr.shape[0]):
			summation[j,:] += arr[frame,j,:]
		summation[j,:] /= arr.shape[0]
	return summation


@njit()
def equidist(p1, p2, npoints=20):
	""" Creates an array of equidistant points between two N-dim points."""
	temp = np.zeros((npoints,p1.shape[0]),dtype=np.float64)
	# loop over x,y,z dimensions
	for i in range(p1.shape[0]):
		temp[:,i] = np.linspace(p1[i],p2[i],npoints)
	return temp

# @njit()
def highsymm_path(symm_points,l,K_step):
	""" Generates high symmetry path 
	along the given points."""

	gamma = np.array([0, 0, 0])
	diff_symm_points = np.diff(symm_points, axis=0)
	path = np.array(symm_points[0],ndmin=2)
	pgp=np.array([0]) if (symm_points[0]==gamma).all() else np.array([],ndmin=1)
	for ii in range(diff_symm_points.shape[0]):
		symmetry_point_linear_displacement = np.max(np.abs(diff_symm_points[ii]))
		steps=int(np.round(symmetry_point_linear_displacement/K_step))
		if (symm_points[ii+1] == gamma).all(): pgp=np.append(pgp,path.shape[0]+steps-1)
		for jj in range(steps):
			path = np.append(path,[path[-1] + diff_symm_points[ii] * 1.0 / steps], axis=0)
	return path, pgp


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


@njit(parallel=True)
def exponential_term(traj,pt):
	traj_0 = mean(traj)
	exponentials = np.zeros((nuq, Natoms), dtype=np.complex128)
	for ii in prange(nuq):
		exponentials[ii] = np.exp(-1j * np.sum(pt[ii] * traj_0, axis=1))
	print("exponentials calculated")
	return exponentials


@njit(parallel=True)
def FT(POS,exponentials,pt):
	'''
	Calculate Fourier transform of position. ONLY for 1 FRAME
	:param traj: in this dimension: POS[natom,3]
	:param pt: high symmetry path
	:return:
	'''
	# print "POS",POS
	kir = np.zeros((nuq,3),dtype=np.complex128)
	if len(POS.shape) != 2:
		raise ValueError

	for ii in prange(nuq):
		kir[ii,0] = Natoms_root_rev*np.sum( POS[:,0] * exponentials[ii])
		kir[ii,1] = Natoms_root_rev*np.sum( POS[:,1] * exponentials[ii])
		kir[ii,2] = Natoms_root_rev*np.sum( POS[:,2] * exponentials[ii])
	return kir


@njit(parallel=True)
def greens_func(traj,pt):
	"""	Takes the Fourier transform of the absolute positions for a given vector.
	Averages all frames and calculates the FT Green's function coeffs at each 
	wave vector q."""
	# R_ka = np.mean(traj,axis=0) # average over all frames

	G_ft = np.zeros((nuq,3,3),dtype=np.complex128) # ka, k'b

	# Calculate exponential term necessarily for FT calculation
	exponentials = exponential_term(traj, pt)
	# For first term
	for fram in range(Nframes):
		Rq      =  FT(traj[fram],exponentials,pt)
		Rq_star =  np.conj(Rq)
		for qq in prange(nuq):
			for alpha in range(3):
				for beta in range(3):
					G_ft[qq,alpha,beta] += Rq[qq,alpha]*Rq_star[qq,beta]
					# print G_ft[qq]
	G_ft=G_ft*(1.0/Nframes)
	# print("Green function first term is done!","Time (seconds): ",time.time()-t_start)
	# For Second term
	# R_mean = np.mean(traj,axis=0)
	# R_mean = np_mean(traj,0)
	R_mean = mean(traj)


	R_mean_q = FT(R_mean,exponentials,pt)
	R_mean_q_star = np.conj(R_mean_q)

	for qq in prange(nuq):
		for alpha in range(3):
			for beta in range(3):
				G_ft[qq,alpha, beta] -= R_mean_q[qq, alpha] * R_mean_q_star[qq, beta]

	# print("Green function Second term is done!","Time (seconds): ",time.time()-t_start)
	print("Green function constructed!")
	print("G_ft.shape=",G_ft.shape)
	return G_ft

# @njit(parallel=True)
def check_hermiticity(G):
	"""	check if G is hermitian 
	!! PROBLEM should check for all atoms separately"""
	for qq in prange(nuq):
		condition_real = (np.round(np.real(G[qq]),5)==np.round(np.real(np.conj(G[qq].T)),5)).all() # redundant
		condition_imag = (np.round(np.imag(G[qq]),5)==np.round(np.imag(np.conj(G[qq].T)),5)).all()
		if (condition_real and condition_imag): # check if G is hermitian
			# print(G[qq])
			print("Matrix is Hermitian, and Determinant is=",np.linalg.det(G[qq]))
		else:
			# print("Matrix is NOT Hermitian\n",np.conj(G)==G)
			print("Matrix is NOT Hermitian for q_n=",qq)
			# print "G.conj is :\n",np.conj(G[qq])
			print("G is :\n",G[qq][0])
			print("G is :\n",G[qq][1])
			print("G is :\n",G[qq][2])
			raise ValueError("Matrix is NOT Hermitian")
			#exit()
	return None

@njit(parallel=True)
def force_constants(G):
	""" Calculates force constants $\Phi_{lk\alpha,l'k'\beta}$ """
	# phi = np.zeros(np.shape(G))

	# check_hermiticity(G)

	Phi = np.zeros(G.shape,dtype=np.complex128) # ka, k'b

	for qq in prange(nuq):
		# print(G[qq])
		Phi[qq] = np.linalg.inv(G[qq])
	####
	return Phi

@njit(parallel=True)
def eigenfreqs(phi_ft,nuq):
	# D = 1/np.sqrt(M*M)* phi_ft
	D = phi_ft
	omega_sq = np.zeros((nuq,3),dtype=np.float64)
	eigenvals_real = np.zeros(3,dtype=np.float64)
	for qq in prange(nuq):
		# eigenvals,eigenvecs = np.linalg.eigh(D[qq])
		eigenvals = np.linalg.eigvals(D[qq])
		eigenvals_real = np.real(eigenvals)
		eidx = eigenvals_real.argsort()[::-1]   # sorting from smallest to largest
		eigenvals_real = eigenvals_real[eidx]
		print("== EIGENVALUES ==\n",eigenvals)
		omega_sq[qq] = eigenvals_real
	print("Success!")
	## Convert to SI units ==>>Hz
	omega_sq *= kB_true*T/(mass*sigma_true*sigma_true)
	print("Frequencies converted to Hz")
	return np.sqrt(omega_sq)

@njit(parallel=False)
def ASR(phi,pgp,nucell):
	"""
	:param phi: is Force matrix at q=0
	:param nucell: is number of atoms in unitcell, For the time being it works only for nucell==1
	:return: ASRed phi_0
	"""
	if nucell != 1: raise ReferenceError

	for ii in pgp:
		phi[ii] = np.imag(phi[ii])*1j  # Zeroing phi at Gamma point (call me Neven, if you complain :)
	return phi


def main():
	t_start = time.time()

	global Natoms
	global Natoms_root_rev
	global Nframes
	global nuq


	load_previous_calculation = False
	load_loaded_traj=False
	# load_previous_calculation = True
	# load_loaded_traj=True

	# set some initial values
	a = np.power(2,(2./3)) # cubic constant in sigma units
	skip_portion =  10 #skip this percent of total time step at the begining
	l = np.array([4, 4, 4])  # lattice size in each direction. THEY MUST BE EQUAL!
	K_step = 2 * np.pi / (a * l[0])
	# print K_step
	# Defining high symmetry points in Kx,Ky,Kz direction ref of the path: http://lampx.tugraz.at/~hadley/ss1/bzones/fcc.php
	gamma = np.array([0,0,0])
	X = np.array([0,2*np.pi/a,0])
	W = np.array([np.pi/a,2*np.pi/a,0])
	K =  np.array([3*np.pi/(2*a),3*np.pi/(2*a),0])
	L = np.array([np.pi/a,np.pi/a,np.pi/a])	
	U = np.array([np.pi/(2*a),2*np.pi/a,np.pi/(2*a)])

	# pt = highsymm_path(np.array([K,gamma,L,W,X,U,X,gamma]),l) # make a path of all points
	# pgp = np.array([1,7])*l[0]  # position of gamma points for ASR  Manually for the time being! ==> Just put index of where gamma points is in pt.
	# plot_ticks = ['K', r'$\Gamma$', 'L', 'W', 'X', 'U', 'X', r'$\Gamma$']
	symm_points = np.array([gamma, X, W, K, gamma, L])
	pt,pgp = highsymm_path(symm_points,l,K_step) # make a path of all points
	# pgp = np.array([0,4])  # position of gamma points for ASR  Manually for the time being! ==> Just put index of where gamma points is in pt.
	plot_ticks = [r'$\Gamma$', 'X', 'W', 'K',r'$\Gamma$', 'L']
	print("pgp=",pgp)

	nuq = pt.shape[0]  # Total number of all points
	print("number of q points=",nuq)


	# !!! new traj system
	traj = np.load(folder_path+trajectory_file.split(".")[0]+'.npy')
	Natoms = traj.shape[1]
	Nframes = traj.shape[0]
	# !!! end new traj system

	Natoms_root_rev = 1.0/np.sqrt(Natoms)
	traj =  traj[int(traj.shape[0]*skip_portion/100):,:,:]
	Nframes = traj.shape[0]

	# MAIN engine
	G_ft = greens_func(traj,pt)     # Calculates green function

	phi_ft = force_constants(G_ft)   # Calculates force matrix in reciprocal space
	phi_ft=ASR(phi_ft,pgp,nucell=1)            # Apply ASR
	freqs = eigenfreqs(phi_ft,nuq)   # Calculates eigen values which is frequencies
	print(" == FREQUENCIES (omega(q)) ==\n",freqs)

	# Save everything, if you wanted to change a little thing in plots
	np.save(folder_path+'temp_pt'+save_flag, pt)
	np.save(folder_path+'temp_freqs'+save_flag,freqs)


	t_end = time.time()
	print("Time taken: {} seconds\n".format(t_end-t_start))

	## project the path to 2D plot
	pt_diff =np.linalg.norm(np.diff(pt,axis=0),axis=1)
	diff_symm_points = np.diff(symm_points, axis=0)
	# print pt_diff
	X=[0]
	plot_ticks_pos = [0]
	for ii in range(pt_diff.shape[0]):
		x=X[ii]+pt_diff[ii]
		X.append(x)
	for ii in range(diff_symm_points.shape[0]):
		plot_ticks_pos.append(plot_ticks_pos[-1]+np.linalg.norm(diff_symm_points[ii]))


	print("High Symmetry Points (projected pos): ",plot_ticks_pos)

	np.savetxt('dispersion.dat',np.array([X, freqs[:, 0]*1e-12,freqs[:, 1]*1e-12,freqs[:, 2]*1e-12]).T)
	
	plt.plot(X, freqs[:, 0]*1e-12,'o-')
	plt.plot(X, freqs[:, 1]*1e-12,'o-')
	plt.plot(X, freqs[:, 2]*1e-12,'o-')
	# print(plot_ticks_pos,plot_ticks)
	plt.xticks(plot_ticks_pos,plot_ticks)
	plt.ylabel('THz')
	plt.savefig(folder_path+"test"+save_flag+".png")
	plt.show()





	return None

if __name__ == '__main__':
	main()