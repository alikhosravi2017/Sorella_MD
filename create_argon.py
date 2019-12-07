#!/usr/bin/env python3
from ase.build import bulk


# latt const https://aip.scitation.org/doi/10.1063/1.1726009
# Ar = bulk('Ar', 'fcc',5.3118)
# =======
unitcell = bulk('Ar', 'fcc',5.3118,orthorhombic=True)
from ase import visualize, Atoms

# Ars= Atoms(bulk, size=(10,10,10))
from ase.build import stack

atoms_a = unitcell
for i in range(10):
    atoms_a = stack(atoms_a, unitcell, axis=0)
atoms_b = atoms_a
for j in range(8):
    atoms_b = stack(atoms_b, atoms_a, axis=1)
all_atoms = atoms_b
for k in range(8):
    all_atoms = stack(all_atoms, atoms_b, axis=2)



visualize.view(all_atoms)
