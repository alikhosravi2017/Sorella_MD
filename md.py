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
np.random.seed(4) # debug


# hardcoded until numba improves support for classes
def initial_parameteres():
    global Natoms,X0, n_a,n_b,n_c,box_sizes,box_half_sizes,Nsteps,cutoff,dump_step,log_step,velocity_zeroing_step,temp_ref,temp_step,Potential_formula,kB_true
    global trajectory_file,trajectory_file_unwrapped, log_file

    Natoms = 0 # will be set, also as global value in premain
    X0=0
    n_a = 3  #Natoms of primitive cell in a direction
    n_b = 3  #Natoms of primitive cell in b direction
    n_c = 3  #Natoms of primitive cell in c direction
    box_sizes = [0,0,0] # from 0 to L in each direction. will be set, also as global value in premain
    box_half_sizes = 0 # will be set, also as global value in premain
    Nsteps = 10**5
    cutoff = 3   #cut off
    dump_step = 100
    log_step = 10
    velocity_zeroing_step =100
    # dt = 0.0005
    temp_ref = 20 # reference tempreature in Kelvin
    temp_step = 100 # thermostat every N steps

    Potential_formula= 'LJ' # Potential type 'LJ' or 'Morse'

    kB_true = 1.38064852e-23  #m2 kg s-2 K-1

    trajectory_file = "traj.xyz"
    trajectory_file_unwrapped = "traj_unwrapped.xyz"
    log_file = "output.dat"

    # constants
    # kBTrue = 1.38064852e-23  #m2 kg s-2 K-1
    # epsilon_True = 1.65e-21 #J
    # sigma_True = 3.4e-10    #m

     
    global force,potential_energy
    global D0,r0,alpha,electronvolt_to_joules,epsilon_true,sigma_true
    global mass,tau,dt
    if Potential_formula == 'Morse':
        D0 =  0.3429 # ev  ## D in Morse potential ## FOR Copper ##
        r0 =  2.866e-10 # meter
        alpha = 1.3588e10 # 1/meter
        # map to LJ units in order to avoid changing the code :)
        electronvolt_to_joules = 1.60218e-19
        epsilon_true = D0 * electronvolt_to_joules  #J
        sigma_true = r0   #meter
        alphar0 = alpha*r0
        force = force_Morse
        potential_energy = potential_energy_Morse
    if Potential_formula == 'LJ':
        epsilon_true = 1.65e-21  #1.977e-21 #J
        sigma_true = 3.4e-10 #3.348e-10    #m
        mass = 6.6335209e-26  # kg
        tau = (1.0)/np.sqrt(epsilon_true/(mass*sigma_true*sigma_true)) # time unit in second
        dt = 0.001 * 1.0/tau * 1e-12
        print("time unit is= {0} Seconds\n and timestep={1} femtoseconds".format(tau,tau*dt*1e15))
        force = force_LJ
        potential_energy = potential_energy_LJ
    return None


@njit(parallel=False)
def pbc(X):
    global box_sizes
    X[:,0] -= box_sizes[0]*np.floor(X[:,0]/box_sizes[0])
    X[:,1] -= box_sizes[1]*np.floor(X[:,1]/box_sizes[1])
    X[:,2] -= box_sizes[2]*np.floor(X[:,2]/box_sizes[2])
    return X

@njit(parallel=True)
def force_LJ(X,F):
    F[:,:] = 0.0
    # numpy vectorization slows down this loOP, CAREFUL!
    for i in prange(Natoms):
        for j in prange(i+1,Natoms):
            # vectorize this
            delta_x = X[i,0]-X[j,0]
            delta_y = X[i,1]-X[j,1]
            delta_z = X[i,2]-X[j,2]
            delta_x += (-1 * box_sizes[0])*np.trunc(delta_x/box_half_sizes[0])
            delta_y += (-1 * box_sizes[1])*np.trunc(delta_y/box_half_sizes[1])
            delta_z += (-1 * box_sizes[2])*np.trunc(delta_z/box_half_sizes[2])
            r2 = delta_x**2 + delta_y**2 + delta_z**2
            if np.sqrt(r2)<cutoff:
                f0= 48*(r2**-7 - 0.5*r2**-4)
                fx = delta_x * f0
                fy = delta_y * f0
                fz = delta_z * f0
                F[i, 0] += fx
                F[i, 1] += fy
                F[i, 2] += fz
                F[j, 0] += -fx
                F[j, 1] += -fy
                F[j, 2] += -fz
    return F

@njit(parallel=True)
def force_Morse(X, F):
    global Natoms,alphar0, box_sizes,box_half_sizes,cutoff
    F[:, :] = 0.0
    for i in prange(Natoms):
        for j in prange(i + 1, Natoms):
            # vectorize this
            delta_x = X[i, 0] - X[j, 0]
            delta_y = X[i, 1] - X[j, 1]
            delta_z = X[i, 2] - X[j, 2]
            delta_x += (-1 * box_sizes[0]) * np.trunc(delta_x / box_half_sizes[0])
            delta_y += (-1 * box_sizes[1]) * np.trunc(delta_y / box_half_sizes[1])
            delta_z += (-1 * box_sizes[2]) * np.trunc(delta_z / box_half_sizes[2])

            r = np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)
            if r < cutoff:
                f0 = -2*alphar0 * (np.exp(-2*alphar0*(r-1)) - 2*alphar0*np.exp(-alphar0*(r-1))) 
                fx = delta_x * f0
                fy = delta_y * f0
                fz = delta_z * f0
                F[i, 0] += fx
                F[i, 1] += fy
                F[i, 2] += fz

                F[j, 0] += -fx
                F[j, 1] += -fy
                F[j, 2] += -fz
    return F



@njit(parallel=False)
def potential_energy_LJ(X):
    global Natoms,box_sizes,box_half_sizes,cutoff
    E = 0.
    for i in prange(Natoms):
        for j in range(i+1,Natoms):
            delta_x = X[i,0]-X[j,0]
            delta_y = X[i,1]-X[j,1]
            delta_z = X[i,2]-X[j,2]
            delta_x += (-1 * box_sizes[0])*np.trunc(delta_x/box_half_sizes[0])
            delta_y += (-1 * box_sizes[1])*np.trunc(delta_y/box_half_sizes[1])
            delta_z += (-1 * box_sizes[2])*np.trunc(delta_z/box_half_sizes[2])

            r = np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
            if r<cutoff:
                E += 4 * ( r**-12 - r**-6)
    return E

@njit(parallel=False)
def potential_energy_Morse(X):
    global Natoms,alphar0, box_sizes,box_half_sizes,cutoff
    E = 0.
    for i in prange(Natoms):
        for j in range(i + 1, Natoms):
            delta_x = X[i, 0] - X[j, 0]
            delta_y = X[i, 1] - X[j, 1]
            delta_z = X[i, 2] - X[j, 2]
            delta_x += (-1 * box_sizes[0]) * np.trunc(delta_x / box_half_sizes[0])
            delta_y += (-1 * box_sizes[1]) * np.trunc(delta_y / box_half_sizes[1])
            delta_z += (-1 * box_sizes[2]) * np.trunc(delta_z / box_half_sizes[2])

            r = np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)
            if r < cutoff:
                E += np.exp(-2*alphar0 * (r - 1)) - 2 * np.exp(-alphar0*(r - 1))
    return E

@njit()
def kinetic_energy(V):
    E = .5*np.sum(V**2)
    return E

@njit()
def temperature(V):
    global Natoms
    return kinetic_energy(V) * 2/(3*Natoms)

@njit()
def thermostat_velocity_rescaling(V):
    global epsilon_true,kB_true,temp_ref
    temp_true = epsilon_true/kB_true # converts to K
    temp_now = temperature(V)*temp_true
    lambda_ = np.sqrt(temp_ref/temp_now)
    V *= lambda_
    return V

@njit(parallel=False)
def velocity_verlet(V,X,F_0):
    global dt,force
    X += dt*V + F_0*dt**2/2
    X = pbc(X) # apply periodic boundary conditions
    F_1 = force(X, F_0)
    V +=  (F_1+F_0) * dt / 2 # calculate velocities at halfstep

    return V,X,F_1

def dump_xyz(f,f3,X,step):
    global Natoms,dump_step

    if step%dump_step!=0:
        return None

    atom_types = Natoms*[1]
    xyz = np.array([atom_types,X[:,0],X[:,1],X[:,2]]).T

    f.write(str(Natoms).encode())
    f.write(b"\n atoms\n")
    np.savetxt(f, xyz,fmt=('%i','%.8f','%.8f','%.8f'))

    dump_xyz_unwrapped(f3,X, step)
    return None


def dump_xyz_unwrapped(f3,X,step):
    global Natoms,box_sizes,box_half_sizes,X0
    X_unwrapped = np.copy(X)
    displacement_M = X-X0
    X_unwrapped[:,0] += (-1 * box_sizes[0]) * np.trunc(displacement_M[:,0] / box_half_sizes[0])
    X_unwrapped[:,1] += (-1 * box_sizes[1]) * np.trunc(displacement_M[:,1] / box_half_sizes[1])
    X_unwrapped[:,2] += (-1 * box_sizes[2]) * np.trunc(displacement_M[:,2] / box_half_sizes[2])

    atom_types = Natoms*[1]
    xyz = np.array([atom_types,X_unwrapped[:,0],X_unwrapped[:,1],X_unwrapped[:,2]]).T

    f3.write(str(Natoms).encode())
    f3.write(b"\n atoms\n")
    np.savetxt(f3, xyz,fmt=('%i','%.8f','%.8f','%.8f'))
    return None


def log(f,X,V,step):
    global epsilon_true,kB_true,log_output,log_step
    global potential_energy
    # print(potential_energy(X))
    if step%log_step!=0:
        return None
    temp_true = epsilon_true/kB_true
    E_kin = kinetic_energy(V)*epsilon_true
    E_pot = potential_energy(X)*epsilon_true
    E_tot = E_kin+E_pot
    temp_now = temperature(V)*temp_true
    log_output = np.array([step,E_kin,E_pot,E_tot,temp_now])
    #print(step,temp_now)
    f.write("\t".join([str(a) for a in log_output])+"\n")
    # np.savetxt(f2, log_output,fmt=('%i','%.8f','%.8f','%.8f','%.8f'))
    return None


def create_atoms(n_a, n_b, n_c):
    unitcell = bulk('Ar', 'fcc', np.power(2,(2./3)),  cubic=True) #,orthorhombic=True)  # ===>  1.122 sigma (from geometry)
    # unitcell = bulk('Ar', 'fcc', 5.4, orthorhombic=True)
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
    print("The cell is= ", all_atoms.get_cell())
    # visualize.view(all_atoms)
    return all_atoms

@njit(parallel=True)
def fix_COM_velocity(V):
    global Natoms
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

def initialize_VXF():
    global n_a,n_b,n_c
    global X0, Natoms
    global box_sizes,box_half_sizes
    atoms = create_atoms(n_a, n_b, n_c)
    Natoms = atoms.get_number_of_atoms()
    print('total number of atoms=', Natoms)
    box_sizes = np.array([atoms.get_cell()[0][0],   atoms.get_cell()[1][1],  atoms.get_cell()[2][2] ] )
    box_half_sizes = box_sizes / 2
    # initialize positions, velocities and forces
    X = atoms.get_positions()
    X0 = np.copy(X)

    # X = np.zeros((Natoms,3))
    V = np.zeros((Natoms,3))
    F = np.zeros((Natoms,3))

    V = np.random.randn(np.shape(V)[0],np.shape(V)[1])
    # print(V)
    V = fix_COM_velocity(V)
    V = thermostat_velocity_rescaling(V)

    return V,X,F

def main():
    global trajectory_file,trajectory_file_unwrapped,log_file
    global Nsteps,temp_step,velocity_zeroing_step
    global force
    global initial_parameteres
    initial_parameteres()

    V,X,F = initialize_VXF()

    #clean old files
    try:
        os.remove(trajectory_file)
        os.remove(trajectory_file_unwrapped)
        os.remove(log_file)
    except OSError:
        pass
    f = open(trajectory_file, "ab")
    f2 = open(trajectory_file_unwrapped, "ab")
    f3 = open(log_file, "a")

    dump_xyz(f,f2,X,0)

    F = force(X,F)

    t_start = time.time()

    for step in range(Nsteps):
        # print(V)
        V,X,F = velocity_verlet(V,X,F)
        dump_xyz(f,f2,X,step)
        log(f3,X,V,step)
        if step % velocity_zeroing_step == 0:
            V = fix_COM_velocity(V)
        if step % temp_step == 0:
            V = thermostat_velocity_rescaling(V)
            #print('yes')
    t_end = time.time()
    print("Time taken: {:.2f} seconds\n".format(t_end-t_start))

    # close files
    f.close()
    f2.close()
    f3.close()
    return None

# object oriented approach, unforunately NUMBA doesnt support classes well so no OOP for now
# class MD(object):
#   def __init__(self,Natoms=20,Nsteps=10**4,box_size=8,dt=10**-5,displacement=2.5,cutoff=2.5,dump_step=10**3,log_step=10**2,thermostat_step=100,temp_ref=160,temp_step=100):
#     self.Natoms = Natoms    #Natoms of atoms
#     self.box_size = box_size # unitless 2D box length
#     self.dt = dt # discretization
#   self.Nsteps =  Nsteps
#     self.displacement = displacement  
#     self.cutoff = cutoff   #cut off
#     self.dump_step = dump_step  
#     self.log_step = log_step  
#     self.thermostat_step = thermostat_step
#     self.temp_ref = temp_step
#     self.temp_step = temp_step    





if __name__ == '__main__':
    main()
