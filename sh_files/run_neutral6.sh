export STAN_NUM_THREADS=6
./stan_models/MPI_neutral_cd4 sample num_warmup=500 num_samples=1000 random seed=1244 data file=datafiles/tester/cd4_data300_s6.Rdump output file=output_csv/tester/neutral6.csv

