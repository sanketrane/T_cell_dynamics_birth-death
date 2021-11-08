#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test_sanket
#SBATCH -D /opt/mesh/eigg/sanket/cmdstan/ontogeny_BUCHI
#SBATCH --nodes=3
#SBATCH --ntasks=149

mpirun -vvvv  MAP_ASM_deltavar_cd4 sample data file=datafiles/cd8_data_s149.Rdump
