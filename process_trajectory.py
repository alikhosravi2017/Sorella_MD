#!/usr/bin/env python3
import time
import numpy as np

folder = './with_lammps/'
fname ='traj_lammps_10K_1fs.xyz'
fpath = folder+fname

def read_xyz(pos=fpath):
	""" Reads .xyz file and returns the trajectory as an array of strings,
	and also returns No. atoms, No. time frames"""
	t_start = time.time()
	Natoms,Nframes = 0,0
	counter = 0
	atoms = []
	with open(pos,"r") as f:
		lines = f.readlines()
		Natoms = int(lines[0])
		Nframes = int(len(lines)/(Natoms+2))
		for i in range(Natoms+2):
			if i>1:
				atoms.append([str(ln) for ln in lines[i].split()][0])
				
	print('Atom types (ordered list):',atoms)
	print("Trajectory read!","Time: {:.2f} (seconds): ".format(time.time()-t_start))
	return lines,Natoms,Nframes

def process_xyz(lines,Natoms,Nframes,save=True):
	"""Takes string array and creates a frame-indexed trajectory array."""
	t_start = time.time()
	traj = np.zeros((Nframes,Natoms,3),dtype=np.float64)
	frame = 0
	k = 0
	counter = 0
	traj_built_flag = True
	atom_type_list = False
	for idx,ln in enumerate(lines):
		if (counter==0 or counter==1):
			k = 0 # atom idx
			counter +=1 
		elif counter>1:
			temp = [float(l) for l in ln.split()]
			traj[frame,k,:] = temp[1:]
			# traj[frame,k,:] = np.fromstring(ln,dtype=np.float64)[1:]
			k += 1
			counter += 1
			if counter==(Natoms+2):
				# print Natoms, ln
				counter = 0
				frame += 1
	print("Trajectory processed!","Time: {:.2f} (seconds): ".format(time.time()-t_start))
	if save:
		""" Saves trajectory array as a numpy binary file """
		t_start = time.time()
		np.save(folder+fname.split(".")[0]+'.npy',traj)
		print("Trajectory processed!","Time: {:.2f} (seconds): ".format(time.time()-t_start))
	return traj


if __name__ == '__main__':
	raw_traj,Natoms,Nframes = read_xyz()
	traj = process_xyz(raw_traj,Natoms,Nframes)

