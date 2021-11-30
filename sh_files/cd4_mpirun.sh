#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test_sanket
#SBATCH -D /opt/mesh/eigg/sanket/cmdstan/ontogeny_BUCHI
#SBATCH --exclude=raasay
#SBATCH --nodes=3
#SBATCH --ntasks=128

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done


echo "stan_models/full_chimera/MAP_${modelname} sample num_samples=500 num_warmup=300 data file=datafiles/cd4_data_s128.Rdump output file=output_csv/full_chimera/ARTF_${modelname}.csv";
mpirun -vvvv stan_models/Full_chimera/MAP_${modelname} sample num_samples=500 num_warmup=300 data file=datafiles/cd4_data_s128.Rdump output file=output_csv/Full_chimera/ARTF_${modelname}.csv


