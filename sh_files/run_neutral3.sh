export STAN_NUM_THREADS=3 
./stan_models/MAP_asm_deltavar_cd4 sample num_warmup=500 num_samples=1000 random seed=1244 data file=datafiles/cd4_data147_s3.Rdump output file=output_csv/asm_deltavar_cd4.csv

