## load libraries
library(cmdstanr)
library(parallel)
library(rstan)

# stan inputs -------------------------------
#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv") 

## prediction points
tseq1 <- seq(14, 120, length.out = 200)
tseq2 <- seq(28, 120, length.out = 200)
tseq3 <- seq(42, 100, length.out = 200)
tseq4 <- seq(70, 160, length.out = 200)
tseq5 <- seq(189, 280, length.out = 200)
ts_pred <- seq(1, 280, length.out = 200)

### data 
data_list <- list(
  K = length(unique(naive_df$mouse_id)),   ## animal level assignments
  J = length(unique(naive_df$group_id)),   ## cohort level assignments
  N = nrow(naive_df),                      ## total number of observations
  
  naive_labelled = naive_df$naive_labelled,  ## counts of labeled naive CD8 T cells
  age_anim = naive_df$age_anim,              ## time points
  
  mouse_id = naive_df$mouse_id,            ## unique id for each animal 
  group_id = naive_df$group_id,            ## unique id for cohort of animals that were stamped together
  t0_group = unique(naive_df$t0_group),     ## animal age when first observation of labeled cells were made
  
  numPred = 200,                 
  tseq1 = tseq1,               
  tseq2 = tseq2,               
  tseq3 = tseq3,               
  tseq4 = tseq4,   
  tseq5 = tseq5,
  ts_pred = ts_pred
)


init_list <- function(){
  list(
    lambda0 = exp(rnorm(1, log(0.03), 0.1)),
    count_bar = rnorm(1, 2e5, 1e3),
    mu_N0 = rnorm(1, 5.5, 0.1),
    sigma_N0 = runif(1, 0, 0.5),
    sigma_counts = runif(1, 0.1, 0.3)
  )
}

parameters_to_plot <- c("lambda",
                        "count_bar",
                        "mu_N0",
                        "sigma_N0",
                        "N0",
                        "sigma_counts"
)

other_rvs <- c("log_lik_counts", "y1_mean_pred", "y2_mean_pred", "y3_mean_pred",
               "y4_mean_pred", "y5_mean_pred")

### params  and interesting quantities to track in the model
parameters <- c(parameters_to_plot, other_rvs)

### compile model
ModelName <- "DDM_N0mouse"
model_file <- file.path("stan_models", paste0(ModelName, '.stan'))
mod <- cmdstan_model(model_file)

### sampling params to fit model to data
stan_fit <- mod$sample(data = data_list, init = init_list, 
                  iter_warmup = 500, iter_sampling = 1000,
                  chains = 4, parallel_chains = 4,
                  save_warmup = T, refresh = 300, 
                  adapt_delta = 0.9)


# saving the stan fit object for individual runs as 
stanfit_save <- rstan::read_stan_csv(stan_fit$output_files())
output_filename <- file.path("Fits_save",  paste0(ModelName, ".rds"))
saveRDS(stanfit_save, file = file.path(output_filename))

# loo-ic values
loo_loglik <- loo::extract_log_lik(stanfit_save, parameter_name = "log_lik_counts", merge_chains = TRUE)
loo_ic <- loo::loo(loo_loglik,  save_psis = FALSE, cores = 4)
ploocv <- loo_ic$estimates

#loo::waic(loo_loglik)

write.csv(ploocv, file = file.path("output_fit", paste0('stats_', ModelName, ".csv")))
