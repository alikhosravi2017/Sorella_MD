#!/usr/bin/env python3
from ase.build import bulk


# latt const https://aip.scitation.org/doi/10.1063/1.1726009
# Ar = bulk('Ar', 'fcc',5.3118)
# =======
unitcell = bulk('Ar', 'fcc',5.3118,orthorhombic=True)
from ase import visualize
from ase.build import stack
from ase.io import write
# Ars= Atoms(bulk, size=(10,10,10))


atoms_a = unitcell
for i in range(10):
    atoms_a = stack(atoms_a, unitcell, axis=0)
atoms_b = atoms_a
for j in range(8):
    atoms_b = stack(atoms_b, atoms_a, axis=1)
all_atoms = atoms_b
for k in range(8):
    all_atoms = stack(all_atoms, atoms_b, axis=2)

# all_atoms.set_cell([all_atoms.get_cell()[0][0], all_atoms.get_cell()[1][1], all_atoms.get_cell()[2][2] ])
write("initial_positions.xyz", all_atoms, format='xyz', parallel=True, append=False)

visualize.view(all_atoms)
print("the cell is=",all_atoms.get_cell())