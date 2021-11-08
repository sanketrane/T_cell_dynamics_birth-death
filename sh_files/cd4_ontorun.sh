#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test_sanket
#SBATCH -D /opt/mesh/eigg/sanket/cmdstan/ontogeny_BUCHI
####SBATCH --nodelist=tiree,raasay,taransay
#SBATCH --ntasks=34

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

echo "stan_models/only_onto/MAP_${modelname}_cd4 sample num_samples=700 num_warmup=300 data file=datafiles/cd4_data_s34.Rdump output file=output_csv/only_onto/${modelname}_cd4.csv";
mpirun -vvvv stan_models/only_onto/MAP_${modelname}_cd4 sample num_samples=700 num_warmup=300 data file=datafiles/cd4_data_s34.Rdump output file=output_csv/only_onto/${modelname}_cd4.csv


