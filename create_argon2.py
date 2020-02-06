from ase import visualize
from ase.build import stack
from ase.io import write
from ase.build import bulk


def create_atoms(n_a, n_b, n_c):
	unitcell = bulk('Ar', 'fcc', 1.5, orthorhombic=True)  # ===> epsilon unit
	atoms_a = unitcell
	for i in range(n_a - 1):
		atoms_a = stack(atoms_a, unitcell, axis=0)
	atoms_b = atoms_a
	for j in range(n_b - 1):
		atoms_b = stack(atoms_b, atoms_a, axis=1)
	all_atoms = atoms_b
	for k in range(n_c - 1):
		all_atoms = stack(all_atoms, atoms_b, axis=2)
	return all_atoms


atoms = create_atoms(4,4, 4)
# global Natoms
# global box_sizes
# Natoms = atoms.get_number_of_atoms()
# box_sizes = np.array([atoms.get_cell()[0][0],   atoms.get_cell()[1][1],  atoms.get_cell()[2][2] ] )


write("geometry.POSCAR",atoms,format="vasp")