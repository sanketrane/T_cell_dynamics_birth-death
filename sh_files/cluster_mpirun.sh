#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test_sanket
#SBATCH -D /opt/mesh/eigg/sanket/cmdstan/ontogeny_BUCHI
#SBATCH --nodes=3
#SBATCH --ntasks=145

mpirun -vvvv stan_models/Full_chimera/MAP_asm_deltavar_cd4 sample data file=datafiles/cd4_data_s145.Rdump
