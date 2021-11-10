library(tidyverse)
library(rstan)

#######
Population <-  'cd4'

## importing data to be fitted 
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


## time points for solving the PDE
dat_time <- ontogeny_data$age.at.S1K    
ont_counts <- ontogeny_data$total_counts
ont_ki <- ontogeny_data$total_kiprop
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
numPred <- length(ts_pred_ont)
numOnt <- length(dat_time)
dat_t0 <- rep(1, numOnt)

for (n_shards in c(numOnt)) {
  stan_rdump(c("numOnt",  "ont_ki", "ont_counts",  "ts_pred_ont",
               "numPred", "n_shards", "dat_time", "dat_t0"),
             file = file.path('datafiles', paste0(Population, '_data_s', n_shards,".Rdump")))
} 


ggplot()+
  geom_point(aes(x=dat_time, y=ont_counts))+
  scale_y_log10()

ggplot(ontogeny_data)+
  geom_point(aes(x=age.at.S1K, y=total_kiprop * 100))+
  scale_y_log10(limits = c(0.3, 100))




#######
Population <-  'cd8'
## importing data to be fitted 
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


## time points for solving the PDE
dat_time <- ontogeny_data$age.at.S1K    
ont_counts <- ontogeny_data$total_counts
ont_ki <- ontogeny_data$total_kiprop
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
t0_ont <- rep(1, 300)
numPred <- length(ts_pred_ont)
numOnt <- length(dat_time)
dat_t0 <- rep(1, numOnt)

for (n_shards in c(numOnt)) {
  stan_rdump(c("numOnt",  "ont_ki", "ont_counts",  "ts_pred_ont",
               "numPred", "n_shards", "dat_time", "dat_t0"),
             file = file.path('datafiles', paste0(Population, '_data_s', n_shards,".Rdump")))
} 


#### Model simulations
rstan::expose_stan_functions("stan_models/only_onto/MAP_asm_rhovar_cd8.stan")
#params_cd4 <- c(336198.622,
#            0.044876653,
#            0.001759622,
#            0.007508088)

params <- c(65226.9368,
            0.012380697,
            0.001809321,
            -0.002370383)
theta <- c(0)
x_r <- c(dat_time[1])
math_reduce(params, theta, 90, dat_t0[1])


artf_dat <- data.frame()
artf_dat <- data.frame("x"= math_reduce(params, theta, ts_pred_ont[1], t0_ont[1]))
for (i in 1:length(ts_pred_ont)) {
  artf_dat[, i] = data.frame("x"=math_reduce(params, theta, ts_pred_ont[i], t0_ont[i]))
}
artf_df <- data.frame(t(artf_dat))

y1_mean <- c()
y2_mean <- c()

set.seed(1357)
for (i in 1:numPred){
  y1_mean[i] = artf_df[i, 1] #+ rnorm(1, 2e6, 2e6)
  y2_mean[i] = artf_df[i, 2] #+ rnorm(1, 0.02, 0.02)
}

ggplot()+
  geom_line(aes(x = ts_pred_ont, y = y1_mean)) +
  geom_point(data = ontogeny_data, aes(x = age.at.S1K, y = total_counts), col = "#094696", size=2, alpha = 0.95) +
  scale_color_manual(name=NULL, values=c(2, 7, 3)) +
  labs(title= paste0('Counts of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5.00, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + guides(col=F)

ggplot()+
  geom_line(aes(x =ts_pred_ont, y = y2_mean)) +
  geom_point(data = ontogeny_data, aes(x = age.at.S1K, y = total_kiprop), col = "#094696", alpha = 0.9, size = 2) +
  scale_color_manual(name=NULL, values=c("#01848a","#094696")) +
  labs(title= paste0('% of Ki67+ cell in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(3, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(0, 1)) +
  myTheme + theme(legend.position = c(0.85,0.85))



