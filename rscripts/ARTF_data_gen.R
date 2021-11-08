library(tidyverse)
library(rstan)

#######
Population <-  'cd4'

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) %>% 
  arrange(age.at.S1K) 

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34)) %>% 
  arrange(age.at.S1K) 

#unique time points in data for odes solver
unique_times_ont <- ontogeny_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_ont <- ontogeny_data$age.at.S1K                        # timepoints in observations 
solve_time_ont <- unique_times_ont$age.at.S1K
time_index_ont <- purrr::map_dbl(data_time_ont, function(x) which(x == solve_time_ont))    # keeping track of index of time point in relation to solve_time


unique_times_chi <- chimera_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_chi <- chimera_data$age.at.S1K 
solve_time_chi <- unique_times_chi$age.at.S1K
time_index_chi <- purrr::map_dbl(data_time_chi, function(x) which(x == solve_time_chi))    # keeping track of index of time point in relation to solve_time


data_time <- c(data_time_ont, data_time_chi)
solve_time <- c(solve_time_ont, solve_time_chi)
dat_t0 <- c(rep(1, length(solve_time_ont)), tb_chi)

## time points for solving the PDE
data_time <- data_fit$age.at.S1K                        # timepoints in observations 
solve_time <- c(unique_times$age.at.S1K)              #unique time points to solve odes  
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time


data_fit <- chimera_data %>% 
  select(contains('S1K'), contains('counts'), contains('kiprop')) %>%
  mutate(dataset = rep('chimera', nrow(chimera_data))) %>%
  bind_rows(ontogeny_data) %>% 
  arrange(age.at.S1K) 

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file)
donor_ki_df <- ki_data %>% filter(subpop == "donor_ki")
host_ki_df <- ki_data %>% filter(subpop == "host_ki")


### importing stan model's functions
rstan::expose_stan_functions("stan_models/Full_chimera/MAP_neutral_cd4.stan")
params <- c(5e5, 0.04, 0.004, 0.005, 0.5, 0.5)
theta <- c(0)
x_i <- c(dat_t0[1:34])
x_r <- c(dat_time[1:34])
test_r <- c()
for (i in 1:34) {
  test_r[i] <- math_reduce(params, theta, x_r[i], 1)
}

math_reduce(params, theta, 5, 1)

pars <- c(0.3, 0.05, 0.004)
init_cond <- c(0,0,3e5,2e5)
solve_chi(100, 60, init_cond, params[2:4])
solve_ode_chi(c(70, 100), c(60, 60), init_cond, pars)
solve_ode_ont(c(5, 100), init_cond[3:4], pars)

## time points for solving the PDE
data_time_ont <- ontogeny_data$age.at.S1K                        # timepoints in observations 
data_time_chi <- chimera_data$age.at.S1K 
tb_chi <- chimera_data$age.at.BMT
ont_counts <- ontogeny_data$total_counts
ont_ki <- ontogeny_data$total_kiprop
chi_counts <- chimera_data$total_counts
N_donor_fraction <- chimera_data$Nfd
donor_ki <- donor_ki_df$ki_prop
host_ki <- host_ki_df$ki_prop
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)
tb_pred1 <- rep(54, 300)
tb_pred2 <- rep(71, 300)
tb_pred3 <- rep(97, 300)
numPred <- length(ts_pred_chi1)
numOnt <- length(data_time_ont)
numChi <- length(data_time_chi)
dat_t0 <- c(rep(1, numOnt), tb_chi)
dat_time <- c(data_time_ont, data_time_chi)
numObs <- length(dat_time)



for (n_shards in c(149)) {
  stan_rdump(c("numOnt",  "ont_ki", "ont_counts", 
               "numChi", "chi_counts",  "N_donor_fraction", "donor_ki", "host_ki",
               "ts_pred_ont", "ts_pred_chi1", "ts_pred_chi2", "ts_pred_chi3",
               "tb_pred1", "tb_pred2", "tb_pred3", "numPred", "n_shards", "dat_time", "dat_t0"),
             file = file.path('datafiles', paste0(Population, '_data_s', n_shards,".Rdump")))
} 


ggplot()+
  geom_point(aes(x=data_time, y=counts))+
  scale_y_log10()

ggplot()+
  geom_point(aes(x=data_time, y=ki_prop * 100))+
  geom_point(data = filter(data_fit, total_counts > 3e7),
             aes(x=age.at.S1K, y=total_kiprop * 100), col=2) +
  scale_y_log10(limits = c(0.3, 100))




#######
Population <-  'cd8'

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) 

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34))

data_fit <- chimera_data %>% 
  select(contains('S1K'), contains('counts'), contains('kiprop')) %>%
  mutate(dataset = rep('chimera', nrow(chimera_data))) %>%
  bind_rows(ontogeny_data) %>% arrange(age.at.S1K) 

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file)
donor_ki_df <- ki_data %>% filter(subpop == "donor_ki")
host_ki_df <- ki_data %>% filter(subpop == "host_ki")

## time points for solving the PDE
time_ont <- data_fit$age.at.S1K                        # timepoints in observations 
time_chi <- chimera_data$age.at.S1K 
tb_chi <- chimera_data$age.at.BMT
total_counts <- data_fit$total_counts
total_ki <- data_fit$total_kiprop
N_donor_fraction <- chimera_data$Nfd
donor_ki <- donor_ki_df$ki_prop
host_ki <- host_ki_df$ki_prop
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)
tb_pred1 <- rep(54, 300)
tb_pred2 <- rep(71, 300)
tb_pred3 <- rep(97, 300)
numPred <- length(ts_pred_chi1)
numOnt <- length(time_ont)
numChi <- length(time_chi)

for (n_shards in c(26, 52)) {
  stan_rdump(c("numOnt", "time_ont", "to_counts","total_ki", 
               "numChi", "time_chi", "tb_chi","N_donor_fraction", "donor_ki", "host_ki",
               "ts_pred_ont", "ts_pred_chi1", "ts_pred_chi2", "ts_pred_chi3",
               "numPred", "tb_pred1", "tb_pred2", "tb_pred3",  "n_shards"),
             file = file.path('datafiles', paste0(Population, '_data', numOnt,
                                                  "_s", n_shards,".Rdump")))
} 


ggplot() +
  geom_point(aes(x=data_time, y=counts)) +
  scale_y_log10()

ggplot()+
  geom_point(aes(x=data_time, y=ki_prop * 100))+
  geom_point(data = filter(data_fit, total_counts > 3e7),
             aes(x=age.at.S1K, y=total_kiprop * 100), col=2) +
  scale_y_log10(limits = c(0.3, 100))



rstan::expose_stan_functions("stan_models/neutral_cd4.stan")
params <- c(5e5, 0.04, 0.004, 0.005, 0.5, 0.5)
theta <- c(0)
x_i <- c(dat_t0[100])
x_r <- c(dat_time[100])
math_reduce(params, theta, x_r, x_i)

init_cond <- c(1e5 * 0.3, 1e5 * 0.7)
map_tester(params, theta, x_r, x_i)
solve_unique(data_time, init_cond, params)

ts_pred <- seq(5, 30)
ode_df <- solve_ode(ts_pred, init_cond, params)
stan_pred_df <- data.frame("time" = ts_pred,
                            "y_pred" = matrix(unlist(ode_df), nrow = length(ode_df), byrow = TRUE))%>%
                                                mutate(total_counts = y_pred.1 + y_pred.2,
                                                       ki_prop = y_pred.1/total_counts)
                                              
set.seed(124)
err_counts <- rnorm(length(ts_pred), 1e4, 1e4)
err_ki <- rnorm(length(ts_pred), 0.01, 0.01)

 
ggplot(stan_pred_df)+
  geom_point(aes(x=time, y=total_counts+ err_counts))+
  scale_y_log10()

ggplot(stan_pred_df)+
  geom_point(aes(x=time, y=(ki_prop + err_ki) * 100))+
  scale_y_log10(limits = c(0.3, 100))

data_time <- ts_pred                      # timepoints in observations 
numObs = length(ts_pred)
counts <- stan_pred_df$total_counts + err_counts
ki_prop <- stan_pred_df$ki_prop + err_ki
ts_pred <- ts_pred
numPred <- length(ts_pred)

for (n_shards in c(1, 3, 6, 12, 26)) {
  stan_rdump(c("numObs", "data_time", "counts","ki_prop", 
               "ts_pred", "numPred", "n_shards"),
             file = file.path('datafiles', "tester", paste0(Population, '_data', numObs,
                                                  "_s", n_shards,".Rdump")))
} 
