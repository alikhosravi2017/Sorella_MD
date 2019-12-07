#!/usr/bin/env python3
from ase.build import bulk


# latt const https://aip.scitation.org/doi/10.1063/1.1726009
Ar = bulk('Ar', 'fcc',5.3118,orthorhombic=True)


print(Ar)