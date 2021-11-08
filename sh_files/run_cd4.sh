export STAN_NUM_THREADS=49
./stan_models/MAP_neutral_cd4 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd4_data147_s49.Rdump output file=output_csv/neutralS49_cd4.csv

./stan_models/MAP_lip_cd4 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd4_data147_s49.Rdump output file=output_csv/lipS49_cd4.csv

./stan_models/MAP_ddm_cd4 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd4_data147_s49.Rdump output file=output_csv/ddmS49_cd4.csv

./stan_models/MAP_rtem_cd4 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd4_data147_s49.Rdump output file=output_csv/rtemS49_cd4.csv

./stan_models/MAP_rtemld_cd4 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd4_data147_s49.Rdump output file=output_csv/rtemldS49_cd4.csv

