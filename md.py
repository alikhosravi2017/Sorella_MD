#!/usr/bin/env python3
import numpy as np
import time

#####################################################
#                                                   #
#              Preliminary code                     #
#       -> needs refactoring into OOP paradigm      #
#                                                   #
#####################################################

Length = 8
displacement = 1.5
Number = 20    #number of atoms
Cutoff = 2.5   #cut off
dump_every = 1000   
log_every = 100  
thermostat_every=100
h = 0.00001
desired_Temperature = 160 #K
kBTrue = 1.38064852e-23  #m2 kg s-2 K-1
epsilon_True = 1.65e-21 #J
sigma_True = 3.4e-10    #m

 #  kB         1 
 #  epsilon    1
 #  sigma      1

def gasdev(idum):
	pass

def Thermostat_velocity_rescaling(V):
	pass

def Force(X,F):
	pass

Force = np.zeros((Number,3,2)) # why a tensor?
KEnergy = np.zeros((Number,2))
PEnergy = np.zeros((Number,2))

# what are these?
X = np.zeros((Number,2))
V = np.zeros((Number,2))
F = np.zeros((Number,2))



# whow does this work?
X_cm, Y_cm = 0,0
T = 100000
k = np.sqrt(Number/2)+1  
m = displacement
N = Number
for j in range(Number/k + 1):
	for i in range(k+1):
		if N!=0:
			X[Number-N,0] = i*m
			X[Number-N,1] = (j+1)*m
			X_cm += X[Number-N,0]/Number
			Y_cm += X[Number-N,1]/Number
			N -= 1

# move CM to the center of the box
X[:,0] += Length/2-X_cm
X[:,1] += Length/2-Y_cm

# dump....?

# what is this?
V[:,0] = gasdev(idum)
V[:,1] = gasdev(idum)

# calculate CM velocity
XMassVelocity, YMassVelocity = 0,0
XMassVelocity += V[:,0]
YMassVelocity += V[:,1]

# set CM velocity to zero
V[:][0] -= XMassVelocity/Number
V[:][1] -= YMassVelocity/Number

Thermostat_velocity_rescaling(V)

# calculate a0 ? acceleration?
Force(X,F)

# MAIN

# clock() ???
# t_start = time.time()

class MD(object):
	def __init__(self):
		self.BoxSize = 8 # unitless 2D box length
		self.displacement = 1.5
		self.Natoms = 20    #number of atoms
		self.Cutoff = 2.5   #cut off
		self.dump_every = 1000   
		self.log_every = 100  
		self.thermostat_every=100
		self.dt = 0.00001 # discretization
		self.desired_Temperature = 160 #K
		

		#constants
		kBTrue = 1.38064852e-23  #m2 kg s-2 K-1
		epsilon_True = 1.65e-21 #J
		sigma_True = 3.4e-10    #m


def main():
	pass