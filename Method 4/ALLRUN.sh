#!/bin/bash

module load anaconda/3-2021.11
module load mpich/ge/gcc/64
source activate mfix-22.2.2
build_mfixsolver --batch --dmp
mpirun -np 24 ./mfixsolver -f vibbubbling.mfx NODESI=6 NODESJ=4 NODESK=1