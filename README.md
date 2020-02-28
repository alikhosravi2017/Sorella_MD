# Sorella MD

An implementations of Molecular Dynamics (MD) in Python3 and Green’s function phonon dispersion calculation directly from MD trajectories as part of an assignment for the course *Computer simulation of condensed matter: from molecular dynamics to quantum Monte Carlo and Tensor Networks* taught by [Prof. Sandro Sorella](https://people.sissa.it/~sorella/) held  at [SISSA](https://www.sissa.it), academic year 2019./2020.

## Description
The MD code can run NVE and NVT ensemble simulations with a velocity-Verlet integrator and a simple velocity rescaling thermostat . While the code can relatively easily be extended to treat multiple atom types, in the present (proof-of-concept)  implementation only monoatomic Lennard-Jones and Morse potentials are implemented. Solid argon was selected as the benchmark system due to its simplicity and plethora of available experimental data for validation.

High performance of the code was achieved with [Numba](http://numba.pydata.org)  just-in-time compilation and parallelisation. 

## File organisation
MD simulation code is located in *md.py*
Phonons are calculated in *md_phonons.py*
Articles references can be found in the *References* subfolder.
All other files are helper files, some for creating the argon lattice, others contain a previous version of the MD code implemented in C++.