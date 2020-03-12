#!/usr/bin/env python3
from ase.build import bulk
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--arepeat", dest = "a_repeat", default = "1", help="Unit cell reparation in \'a\' direction")
parser.add_argument("-b", "--brepeat",dest ="b_repeat",  default = "1", help="Unit cell reparation in \'b\' direction")
parser.add_argument("-c", "--crepeat",dest ="c_repeat",  default = "1", help="Unit cell reparation in \'c\' direction")
parser.add_argument("-n", "--nunitcell",dest ="n_unitcell",  default = "1", help="Unit cell reparation in all 3 direction")
args = parser.parse_args()
a_repeat = int(args.a_repeat)
b_repeat = int(args.b_repeat)
c_repeat = int(args.c_repeat)
n_unitcell = int(args.n_unitcell)
if a_repeat==b_repeat==c_repeat==1:
    a_repeat = b_repeat = c_repeat = n_unitcell
# latt const https://aip.scitation.org/doi/10.1063/1.1726009   ==>> angstrom
# Ar = bulk('Ar', 'fcc',5.3118)
# =======
unitcell = bulk('Ar', 'fcc',1.5,orthorhombic=True)  #===> epsilon unit
from ase import visualize
from ase.build import stack
from ase.io import write
import numpy as np
# Ars= Atoms(bulk, size=(10,10,10))


atoms_a = unitcell
for i in range(a_repeat-1):
    atoms_a = stack(atoms_a, unitcell, axis=0)
atoms_b = atoms_a
for j in range(b_repeat-1):
    atoms_b = stack(atoms_b, atoms_a, axis=1)
all_atoms = atoms_b
for k in range(c_repeat-1):
    all_atoms = stack(all_atoms, atoms_b, axis=2)

# write("initial_positions.xyz", all_atoms, format='xyz', parallel=True, append=False)
pos_ = all_atoms.get_positions()

# visualize.view(all_atoms)
# print("the cell is=",all_atoms.get_cell())
cell_ = np.array([[all_atoms.get_cell()[0][0],all_atoms.get_cell()[1][1],all_atoms.get_cell()[2][2]]])
cell_plus_pos = np.append(cell_,pos_,axis=0)
# print cell_
np.savetxt("initial_positions_{0}_atoms.xyz".format(len(pos_)),cell_plus_pos)
