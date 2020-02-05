#!/usr/bin/env python3
import numpy as np
import time
import os
from numba import njit,prange
from ase import visualize
from ase.build import stack
from ase.io import write
import numpy as np
from ase.build import bulk
#####################################################
#                                                   #
#              Preliminary code                     #
#       -> needs refactoring into OOP paradigm      #
#                                                   #
#####################################################
Natoms = 0 # will be set, also as global value in premain
n_a = 2  #Natoms of primitive cell in a direction
n_b = 2  #Natoms of primitive cell in b direction
n_c = 2  #Natoms of primitive cell in c direction
box_sizes = [0,0,0] # from 0 to L in each direction. will be set, also as global value in premain
Nsteps = 10**6
cutoff = 3   #cut off
dump_step = 1000
log_step = 1000
velocity_zeroing_step =100
dt = 0.00001
temp_ref = 1 # reference tempreature in Kelvin
temp_step = 1000 # thermostat every N steps
kB_true = 1.38064852e-23  #m2 kg s-2 K-1
epsilon_true = 1.977e-21 #J
sigma_true = 3.348e-10    #m
# random_seed = 8
trajectory_file = "traj.xyz"
log_file = "output.dat"


np.random.seed(8)

try:
    os.remove(trajectory_file)
    os.remove(log_file)
except OSError: pass
f = open(trajectory_file, "ab")
f2 = open(log_file, "a")
#  kB         1
#  epsilon    1
#  sigma      1



def pbc(X):
    X[:,0] -= box_sizes[0]*np.floor(X[:,0]/box_sizes[0])
    X[:,1] -= box_sizes[1]*np.floor(X[:,1]/box_sizes[1])
    X[:,2] -= box_sizes[2]*np.floor(X[:,2]/box_sizes[2])
    return X

# needs refactoring
@njit(parallel=True)
def force(X,F):
    #print("== X ==\n",X)
    F[:,:] = 0.0

    box_half_sizes = box_sizes/2

    # this loop can be vectorized --see numpy documentation
    for i in prange(Natoms):
        for j in prange(Natoms):
            # vectorize this
            if j!=i:
        #what is this?
                solid_distance_x = np.abs(X[i,0]-X[j,0])
                solid_distance_y = np.abs(X[i,1]-X[j,1])
                solid_distance_z = np.abs(X[i,2]-X[j,2])
                delta_x = -solid_distance_x + box_sizes[0] * np.floor(solid_distance_x/box_half_sizes[0])
                delta_y = -solid_distance_y + box_sizes[1] * np.floor(solid_distance_y/box_half_sizes[1])
                delta_z = -solid_distance_z + box_sizes[2] * np.floor(solid_distance_z/box_half_sizes[2])
                r2 = delta_x**2 + delta_y**2 + delta_z**2

                # what is this
                if np.sqrt(r2)<cutoff:
                    if X[i,0]<X[j,0]:
                        F[i,0] += delta_x*48*((1/r2**7) - 1/(2*r2**4))
                    else:
                        F[i,0] += -delta_x*48*((1/r2**7) - 1/(2*r2**4))

                    if X[i,1]<X[j,1]:
                        F[i,1] += delta_y*48*((1/r2**7) - 1/(2*r2**4))
                    else:
                        F[i,1] += -delta_y*48*((1/r2**7) - 1/(2*r2**4))

                    if X[i,2]<X[j,2]:
                        F[i,2] += delta_x*48*((1/r2**7) - 1/(2*r2**4))
                    else:
                        F[i,2] += -delta_x*48*((1/r2**7) - 1/(2*r2**4))
    # F = -F.T + F + F[np.diag_indices_from(F)]
    return F

# needs to be vectorized
@njit(parallel=False)
def potential_energy(X):
    E = 0
    for i in range(Natoms):
        for j in range(i+1,Natoms):
            if np.abs(X[i,0]-X[j,0]) > box_sizes[0]/2:
                delta_x = box_sizes[0] - np.abs(X[i,0]-X[j,0])
            else:
                delta_x = np.abs(X[i,0]-X[j,0])

            if np.abs(X[i,1]-X[j,1]) > box_sizes[1]/2:
                delta_y = box_sizes[1] - np.abs(X[i,1]-X[j,1])
            else:
                delta_y = np.abs(X[i,1]-X[j,1])
            
            if np.abs(X[i,2]-X[j,2]) > box_sizes[2]/2:
                delta_z = box_sizes[2] - np.abs(X[i,2]-X[j,2])
            else:
                delta_z = np.abs(X[i,2]-X[j,2])

            r = np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
            if r<cutoff:
                E += 4 * ( r**-12 - r**-6)
    return E


def kinetic_energy(V):
    # is this ok?
    E = .5*np.sum(V**2)
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
    xyz = np.array([atom_types,X[:,0],X[:,1],X[:,2]]).T

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
    print(log_output)
    f2.write("\t".join([str(a) for a in log_output])+"\n")
    # np.savetxt(f2, log_output,fmt=('%i','%.8f','%.8f','%.8f','%.8f'))
    return None

#Forces = np.zeros((Natoms,3,2)) # why a tensor?
# KEnergy = np.zeros((Natoms,2))
# PEnergy = np.zeros((Natoms,2))

### positioning
def create_atoms(n_a, n_b, n_c):
    unitcell = bulk('Ar', 'fcc', np.power(2,(1./6))*2/np.sqrt(2), cubic=True)  # ===>  1.122 sigma
    # visualize.view(unitcell)
    # print unitcell.get_cell()
    atoms_a = unitcell
    for i in range(n_a - 1):
        atoms_a = stack(atoms_a, unitcell, axis=0)
    atoms_b = atoms_a
    for j in range(n_b - 1):
        atoms_b = stack(atoms_b, atoms_a, axis=1)
    all_atoms = atoms_b
    for k in range(n_c - 1):
        all_atoms = stack(all_atoms, atoms_b, axis=2)
    # print all_atoms.get_cell()
    return all_atoms
###

@njit(parallel=True)
def fix_COM_velocity(V):
    # calculate CM velocity
    # XMassVelocity, YMassVelocity = 0,0
    XMassVelocity = np.sum(V[:, 0])
    YMassVelocity = np.sum(V[:, 1])
    ZMassVelocity = np.sum(V[:, 2])

    # set CM velocity to zero
    V[:, 0] -= XMassVelocity / Natoms
    V[:, 1] -= YMassVelocity / Natoms
    V[:, 2] -= ZMassVelocity / Natoms
    return  V

# @njit(parallel=True)
def pre_main():
    atoms = create_atoms(n_a, n_b, n_c)
    global Natoms
    global box_sizes
    Natoms = atoms.get_number_of_atoms()
    print('total number of atoms=', Natoms)
    box_sizes = np.array([atoms.get_cell()[0][0],   atoms.get_cell()[1][1],  atoms.get_cell()[2][2] ] )

    # initialize positions, velocities and forces
    X = atoms.get_positions()
    # X = np.zeros((Natoms,3))
    V = np.zeros((Natoms,3))
    F = np.zeros((Natoms,3))

    V = np.random.randn(np.shape(V)[0],np.shape(V)[1])
    V = fix_COM_velocity(V)
    V = thermostat_velocity_rescaling(V)
    return V,X,F

def main():

    V,X,F = pre_main()
    dump_xyz(X,0,trajectory_file)
    # print(XMassVelocity,YMassVelocity)


    # calculate a0 ? acceleration?
    F = force(X,F)

    t_start = time.time()

    for step in range(Nsteps):
        # print(V)
        V,X,F = velocity_verlet(V,X,F)
        dump_xyz(X,step,trajectory_file)
        log(X,V,step,log_file)
        if step % velocity_zeroing_step == 0:
            V = fix_COM_velocity(V)
        if step % temp_step == 0:
            V = thermostat_velocity_rescaling(V)
            #print('yes')
    t_end = time.time()
    #print("Time taken: {:.2f} seconds\n".format(t_end-t_start))

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
