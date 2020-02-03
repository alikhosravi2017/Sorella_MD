#!/usr/bin/env python3
import numpy as np
import time
import os
from numba import njit,prange
#####################################################
#                                                   #
#              Preliminary code                     #
#       -> needs refactoring into OOP paradigm      #
#                                                   #
#####################################################
Natoms = 50    #Natoms of atoms
box_size = 8
Nsteps = 10**6
displacement = 1.5
cutoff = 2.5   #cut off
dump_step = 1000   
log_step = 1000  
dt = 0.00001
temp_ref = 60 #K
temp_step = 100 # ever 100th step
kB_true = 1.38064852e-23  #m2 kg s-2 K-1
epsilon_true = 1.65e-21 #J
sigma_true = 3.4e-10    #m
# random_seed = 8
trajectory_file = "traj.xyz"
log_file = "output.dat"


np.random.seed(8)

# clean old files
os.remove(trajectory_file)
os.remove(log_file)
# open files for writting
f = open(trajectory_file, "ab")
f2 = open(log_file, "a")


 #  kB         1 
 #  epsilon    1
 #  sigma      1



def pbc(X):
	X -= box_size*np.floor(X/box_size)
	# X[:,0] -= box_size*np.floor(X[:,0]/box_size)
	# X[:,1] -= box_size*np.floor(X[:,1]/box_size)
	return X

# needs refactoring
@njit(parallel=True)
def force(X,F):
	F[:,:] = 0

	box_half_size = box_size/2
	# this loop can be vectorized --see numpy documentation
	for i in prange(Natoms):
		for j in prange(Natoms):
			# vectorize this
			if j!=i:
			#what is this?
				solid_distance_x = np.abs(X[i,0]-X[j,0])
				solid_distance_y = np.abs(X[i,1]-X[j,1])
				delta_x = -solid_distance_x + box_size * np.floor(solid_distance_x/box_half_size)
				delta_y = -solid_distance_y + box_size * np.floor(solid_distance_y/box_half_size)
				r2 = delta_x**2 + delta_y**2
		
				# what is this
				if np.sqrt(r2)<cutoff:
					if X[i,0]<X[j,0]:
						F[i,0] += delta_x*48*(1/r2**7) - 1/(2*r2**4)
					else:
						F[i,0] += -delta_x*48*(1/r2**7) - 1/(2*r2**4)
					if X[i,1]<X[j,1]:
						F[i,1] += delta_y*48*(1/r2**7) - 1/(2*r2**4)
					else:
						F[i,1] += -delta_y*48*(1/r2**7) - 1/(2*r2**4)
		
	# F = -F.T + F + F[np.diag_indices_from(F)]

	return F

# needs to be vectorized
@njit(parallel=False)
def potential_energy(X):
	E = 0
	for i in range(Natoms):
		for j in range(i+1,Natoms):
			if np.abs(X[i,0]-X[j,0]) > box_size/2:
				delta_x = box_size - np.abs(X[i,0]-X[j,0])		
			else:
				delta_x = np.abs(X[i,0]-X[j,0])

			if np.abs(X[i,1]-X[j,1]) > box_size/2:
				delta_y = box_size - np.abs(X[i,1]-X[j,1])
			else:
				delta_y = np.abs(X[i,1]-X[j,1])
			r = np.sqrt(delta_x**2 + delta_y**2)
			if r<cutoff:
				E += 4 * ( r**-12 - r**-6)
	# for i in range(Natoms):
	# 	for j in range(Natoms):
	# 		if j!=i:
	# 			if np.abs(X[i,0]-X[j,0]) > box_size/2:
	# 				delta_x = box_size - np.abs(X[i,0]-X[j,0])		
	# 			else:
	# 				delta_x = np.abs(X[i,0]-X[j,0])
	
	# 			if np.abs(X[i,1]-X[j,1]) > box_size/2:
	# 				delta_y = box_size - np.abs(X[i,1]-X[j,1])
	# 			else:
	# 				delta_y = np.abs(X[i,1]-X[j,1])
	# 			r = np.sqrt(delta_x**2 + delta_y**2)
	# 			if r<cutoff:
	# 				E += 4 * ( r**-12 - r**-6)
# removed bcs Uji =0
	# E *= .5 # because Uij=Uji
	return E


def kinetic_energy(V):
	#E = np.zeros((Natoms,2))
	E = 0
	# is this ok?
	E = .5*(np.sum(V[:,0]**2) + np.sum(V[:,1]**2))
	return E


def temperature(V):
	#is this ok
    return kinetic_energy(V) * 2/(3*Natoms)


def thermostat_velocity_rescaling(V):
	temp_true = epsilon_true/kB_true # converts to K
	temp_now = temperature(V)*temp_true
	lambda_ = np.sqrt(temp_ref/temp_now)
	V *= lambda_
	return V


def velocity_verlet(V,X,F_0):
	#verlet loop
	# V[:,0] += dt/2*F[:,0]
	# V[:,1] += dt/2*F[:,1]
	# X[:,0] += dt*V[:,0]
	# X[:,1] += dt*V[:,1]
	V += dt/2 * F_0
	X += dt*V
	X = pbc(X)

	# calculate accelerations
	F_1 = np.array(force(X,F_0))
	# calculate velocities at halfstep
	# V[:,0] += dt/2 * F[:,0]
	# V[:,1] += dt/2 * F[:,1]
	V += dt/2*F_1
	F_0 = F_1
	return V,X,F_0

def dump_xyz(X,step,fname):
	if step%dump_step!=0:
		return None

	atom_types = Natoms*[1]
	xyz = np.array([atom_types,X[:,0],X[:,1],np.zeros(np.shape(X)[0])]).T

	f.write(str(Natoms).encode())
	f.write(b"\n atoms\n")
	np.savetxt(f, xyz,fmt=('%i','%.8f','%.8f','%.8f'))
	return None

def log(X,V,step,fname):
	if step%log_step!=0:
		return None

	temp_true = epsilon_true/kB_true	
	E_kin = kinetic_energy(V)*epsilon_true
	E_pot = potential_energy(X)*epsilon_true
	E_tot = E_kin+E_pot
	temp_now = temperature(V)*temp_true
	log_output = np.array([step,E_kin,E_pot,E_tot,temp_now])
	# print(log_output)

	f2.write("\t".join([str(a) for a in log_output])+"\n")
	# np.savetxt(f2, log_output,fmt=('%i','%.8f','%.8f','%.8f','%.8f'))
	return None

#Forces = np.zeros((Natoms,3,2)) # why a tensor?
# KEnergy = np.zeros((Natoms,2))
# PEnergy = np.zeros((Natoms,2))

@njit(parallel=True)
def pre_main():


	# initialize positions, velocities and forces
	X = np.zeros((Natoms,2))
	V = np.zeros((Natoms,2))
	F = np.zeros((Natoms,2))



	# how does this work?
	X_cm, Y_cm = 0,0

	k = np.sqrt(Natoms/2)+1

	m = displacement
	N = Natoms

	for j in range(int(Natoms/k)):
		for i in range(int(k)+1):
			# strange loop...
			if N!=0:
				X[Natoms-N,0] = i*m
				X[Natoms-N,1] = (j+1)*m
				X_cm += X[Natoms-N,0]/Natoms
				Y_cm += X[Natoms-N,1]/Natoms
				N -= 1

	# move CM to the center of the box
	X[:,0] += box_size/2-X_cm
	X[:,1] += box_size/2-Y_cm




	V = np.random.randn(np.shape(V)[0],np.shape(V)[1])

	# calculate CM velocity
	# XMassVelocity, YMassVelocity = 0,0
	XMassVelocity = np.sum(V[:,0])
	YMassVelocity = np.sum(V[:,1])

	# set CM velocity to zero
	V[:,0] -= XMassVelocity/Natoms
	V[:,1] -= YMassVelocity/Natoms


	XMassVelocity = np.sum(V[:,0])
	YMassVelocity = np.sum(V[:,1])



	# MAIN

	return V,X,F

def main():

	V,X,F = pre_main()

	dump_xyz(X,0,trajectory_file)

	# print(XMassVelocity,YMassVelocity)
	V = thermostat_velocity_rescaling(V)

	# calculate a0 ? acceleration?
	F = force(X,F)

	t_start = time.time()

	for step in range(Nsteps):
		# print(V)
		V,X,F = velocity_verlet(V,X,F)
		dump_xyz(X,step,trajectory_file)
		log(X,V,step,log_file)
		if step%temp_step==0:
			V = thermostat_velocity_rescaling(V)
	t_end = time.time()
	print("Time taken: {:.2f} seconds\n".format(t_end-t_start))

	# close files
	f.close()
	f2.close()
	return None

# class MD(object):
# 	def __init__(self,Natoms=20,Nsteps=10**4,box_size=8,dt=10**-5,displacement=2.5,cutoff=2.5,dump_step=10**3,log_step=10**2,thermostat_step=100,temp_ref=160,temp_step=100):
# 		self.Natoms = Natoms    #Natoms of atoms
# 		self.box_size = box_size # unitless 2D box length
# 		self.dt = dt # discretization
#		self.Nsteps =  Nsteps
# 		self.displacement = displacement	
# 		self.cutoff = cutoff   #cut off
# 		self.dump_step = dump_step  
# 		self.log_step = log_step  
# 		self.thermostat_step = thermostat_step
# 		self.temp_ref = temp_step
# 		self.temp_step = temp_step		


# 		#constants
# 		kBTrue = 1.38064852e-23  #m2 kg s-2 K-1
# 		epsilon_True = 1.65e-21 #J
# 		sigma_True = 3.4e-10    #m


# def main():
# 	#MD()
# 	pass

if __name__ == '__main__':
	main()