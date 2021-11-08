export STAN_NUM_THREADS=52
./stan_models/MAP_neutral_cd8 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd8_data156_s52.Rdump output file=output_csv/neutralS52_cd8.csv

./stan_models/MAP_lip_cd8 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd8_data156_s52.Rdump output file=output_csv/lipS52_cd8.csv

./stan_models/MAP_ddm_cd8 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd8_data156_s52.Rdump output file=output_csv/ddmS52_cd8.csv

./stan_models/MAP_rtem_cd8 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd8_data156_s52.Rdump output file=output_csv/rtemS52_cd8.csv

./stan_models/MAP_rtemld_cd8 sample num_warmup=300 num_samples=1000 random seed=1244 data file=datafiles/cd8_data156_s52.Rdump output file=output_csv/rtemldS52_cd8.csv


