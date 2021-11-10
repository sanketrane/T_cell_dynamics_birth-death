## load libraries
library(cmdstanr)
library(parallel)
library(rstan)

# stan inputs -------------------------------
#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv")

## prediction points
tseq1 <- seq(14, 120, length.out = 100)
tseq2 <- seq(28, 120, length.out = 100)
tseq3 <- seq(42, 100, length.out = 100)
tseq4 <- seq(70, 160, length.out = 100)
tseq5 <- seq(189, 280, length.out = 100)

### data 
data_list <- list(
  K = length(unique(naive_df$mouse_id)),   ## animal level assignments
  J = length(unique(naive_df$group_id)),   ## cohort level assignments
  N = nrow(naive_df),                      ## total number of observations
  
  naive_labelled = naive_df$naive_labelled,  ## counts of labeled naive CD8 T cells
  age_anim = naive_df$age_anim,              ## time points
  
  mouse_id = naive_df$mouse_id,             ## unique id for each animal 
  group_id = naive_df$group_id,             ## unique id for cohort of animals that were stamped together
  
  day_stamp = unique(naive_df$day_stamp),   ## animal age at time-stamp 
  t0_group = unique(naive_df$t0_group),     ## animal age when first observation of labeled cells were made
  
  numPred = 100,                 
  tseq1 = tseq1,               
  tseq2 = tseq2,               
  tseq3 = tseq3,               
  tseq4 = tseq4,   
  tseq5 = tseq5
)

init_list <- function(){
  list(
    mu_lambda = exp(rnorm(1, log(0.05), 0.1)),
    mu_N0     = rnorm(1, 5.5, 0.1),
    r_lamda   = exp(rnorm(1, log(0.01), 0.1)),
    sigma_counts = runif(1, 0.1, 0.3),
    sigma_N0   = runif(1, 0, 0.5),
    sigma_lambda    = runif(1, 0, 0.01)
  )
}

parameters_to_plot <- c("r_lambda",
                        "N0",
                        "mu_N0",
                        "sigma_N0",
                        "lambda",
                        "mu_lambda",
                        "sigma_lambda",
                        "sigma_counts"
)

other_rvs <- c("log_lik_counts", "y1_mean_pred", "y2_mean_pred", "y3_mean_pred",
               "y4_mean_pred", "y5_mean_pred")

### params  and interesting quantities to track in the model
parameters <- c(parameters_to_plot, other_rvs)

### compile model
ModelName <- "PureASM_N0grp_L0mouse"
model_file <- file.path("stan_models", paste0(ModelName, '.stan'))
mod <- cmdstan_model(model_file)

### sampling params to fit model to data
stan_fit <- mod$sample(data = data_list, init = init_list, 
                  iter_warmup = 500, iter_sampling = 1000,
                  chains = 4, parallel_chains = 4,
                  save_warmup = T, refresh = 100, 
                  adapt_delta = 0.9)

# saving the stan fit object for individual runs as 
stanfit_save <- rstan::read_stan_csv(stan_fit$output_files())
output_filename <- file.path("Fits_save",  paste0(ModelName, ".rds"))
write_rds(stanfit_save, file = file.path(output_filename))


### parameters table
num_pars <- 6 + data_list$K + data_list$J

ptable <- monitor(as.array(stanfit_save, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
write.csv(out_table, file = file.path("output_fit", paste0('params_', ModelName, ".csv")))


# loo-ic values
loo_loglik <- loo::extract_log_lik(stanfit_save, parameter_name = "log_lik_counts", merge_chains = TRUE)
loo_ic <- loo::loo(loo_loglik,  save_psis = FALSE, cores = 4)
ploocv <- loo_ic$estimates

#loo::waic(loo_loglik)

write.csv(ploocv, file = file.path("output_fit", paste0('stats_', ModelName, ".csv")))
