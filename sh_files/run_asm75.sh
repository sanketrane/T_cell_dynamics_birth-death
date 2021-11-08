export STAN_NUM_THREADS=70
./stan_models/MAP_asm_deltavar_cd8 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd8_data156_s78.Rdump output file=output_csv/asm_deltavar_cd8.csv

