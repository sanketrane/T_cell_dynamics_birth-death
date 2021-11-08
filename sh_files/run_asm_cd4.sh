export STAN_NUM_THREADS=149
./stan_models/Full_chimera/MAP_asm_rhovar_cd4 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd4_data_s149.Rdump output file=output_csv/asm_FullS149_rhovar_cd4.csv

