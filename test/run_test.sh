#!/bin/bash

##  Check the program output given the test set of parameters in paramteres.dat

mpirun -n 4 ../../optical_flux_lattice_diag --in-path ./ -d 1 -v 3 --out-path ./ -x 4 -y 4 --nbr 4 --method 0 --eigenvalues-file 1

mpirun -n 4 ../../optical_flux_lattice_diag --in-path ./ -d 1 -v 3 --out-path ./ -x 5 -y 5 --nbr 5 --method 1 --eigenvalues-file 1

mpirun -n 4 ../../optical_flux_lattice_diag --in-path ./ -d 1 -v 3 --out-path ./ -x 5 -y 4 --nbr 5 --method 1 --basis 1 --eigenvalues-file 1

mpirun -n 4 ../../single_particle_model --in-path ./  -v 3 --out-path ./ --plot-x 6 --plot-y 6 --x-cut 10 --y-cut 10 --plot-band-width 1 
