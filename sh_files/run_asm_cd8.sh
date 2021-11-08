export STAN_NUM_THREADS=52
./stan_models/MAP_asm_rhovar_cd8 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd8_data156_s52.Rdump output file=output_csv/asmS52_rhovar_cd8.csv

