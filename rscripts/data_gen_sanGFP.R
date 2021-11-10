library(tidyverse)
library(rstan)

#######
Population <-  'cd4'

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_new.csv"))  
chimera_data <- read.csv(chimera_file) %>%
  arrange(age.at.S1K) %>%
  filter(Nfd < 1.0)

unique_times_chi <- chimera_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_chi <- chimera_data$age.at.S1K 
solve_time_chi <- unique_times_chi$age.at.S1K
time_index_chi <- purrr::map_dbl(data_time_chi, function(x) which(x == solve_time_chi))    # keeping track of index of time point in relation to solve_time

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file) %>%
  right_join(chimera_data) %>% 
  select(-contains("total"), -contains("Nfd"))
donor_ki_df <- ki_data %>% filter(subpop == "donor_ki")
host_ki_df <- ki_data %>% filter(subpop == "host_ki")


#### importing data for GFP positive proportions
GFPposKipos_df <- readxl::read_excel(file.path("datafiles", "RagGFP_ontogeny_pooled.xlsx"), sheet = 6)%>%
  select('mouseID', 'age', 'LN.4nai', 'SP.4nai')


GFPki_df <- readxl::read_excel(file.path("datafiles", "RagGFP_ontogeny_pooled.xlsx"), sheet = 5)%>%
  select('mouseID', 'age','LN.4nai', 'SP.4nai') %>%
  left_join(GFPposKipos_df, by = c('mouseID', 'age')) %>%
  mutate(Ki_percent = ((LN.4nai.x + LN.4nai.y) + (SP.4nai.x + SP.4nai.y))/2,
         total_kiprop = Ki_percent/100,
         age.at.S1K = age) %>%
  select("mouseID", contains("S1K"), contains("total"))

GFP_data <- readxl::read_excel(file.path("datafiles", "RagGFP_ontogeny_pooled.xlsx"), sheet = 1)%>%
  select('mouseID', 'age', 'SP.4nai', 'LN.4nai') %>%
  mutate(total_counts = LN.4nai + SP.4nai,
         age.at.S1K = age) %>%
  select("mouseID", contains("S1K"), contains("total")) %>%
  left_join(GFPki_df, by = c('mouseID', 'age.at.S1K'))


ggplot() +
  geom_point(data=chimera_data, aes(x= age.at.S1K, y = total_counts), col=2)+
  geom_point(data=GFP_data, aes(x= age.at.S1K, y = total_counts), col=4) +
  scale_y_log10() + scale_x_log10()

ggplot() +
  geom_point(data=chimera_data, aes(x= age.at.S1K, y = total_kiprop*100), col=2)+
  geom_point(data=GFP_data, aes(x= age.at.S1K, y = total_kiprop*100), col=4) +
  scale_y_log10() + scale_x_log10()

## time points for solving the PDE
data_time_ont <- GFP_data$age.at.S1K                        # timepoints in observations 
data_time_chi <- chimera_data$age.at.S1K 
tb_chi <- chimera_data$age.at.BMT
ont_counts <- GFP_data$total_counts
ont_ki <- GFP_data$total_kiprop
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

for (n_shards in c(numObs)) {
  stan_rdump(c("numOnt",  "ont_ki", "ont_counts", 
               "numChi", "chi_counts",  "N_donor_fraction", "donor_ki", "host_ki",
               "ts_pred_ont", "ts_pred_chi1", "ts_pred_chi2", "ts_pred_chi3",
               "tb_pred1", "tb_pred2", "tb_pred3", "numPred", "n_shards", "dat_time", "dat_t0"),
             file = file.path('datafiles', paste0(Population, '_data_s', n_shards,".Rdump")))
} 

########################################################################################################################################################################
########################################################################################################################################################################

Population <-  'cd8'

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_new.csv"))  
chimera_data <- read.csv(chimera_file) %>%
  arrange(age.at.S1K) %>%
  filter(Nfd < 1.0)

unique_times_chi <- chimera_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_chi <- chimera_data$age.at.S1K 
solve_time_chi <- unique_times_chi$age.at.S1K
time_index_chi <- purrr::map_dbl(data_time_chi, function(x) which(x == solve_time_chi))    # keeping track of index of time point in relation to solve_time

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file) %>%
  right_join(chimera_data) %>% 
  select(-contains("total"), -contains("Nfd"))
donor_ki_df <- ki_data %>% filter(subpop == "donor_ki")
host_ki_df <- ki_data %>% filter(subpop == "host_ki")

#### importing data for GFP positive proportions
GFPposKineg_df <- readxl::read_excel(file.path("datafiles", "RagGFP_ontogeny_pooled.xlsx"), sheet = 3)%>%
  select('mouseID', 'age', 'LN.8nai')
GFPposKipos_df <- readxl::read_excel(file.path("datafiles", "RagGFP_ontogeny_pooled.xlsx"), sheet = 6)%>%
  select('mouseID', 'age', 'LN.8nai')


GFPki_df <- readxl::read_excel(file.path("datafiles", "RagGFP_ontogeny_pooled.xlsx"), sheet = 5)%>%
  select('mouseID', 'age','LN.8nai') %>%
  left_join(GFPposKipos_df, by = c('mouseID', 'age')) %>%
  mutate(Ki_percent =LN.8nai.x + LN.8nai.y,
         total_kiprop = Ki_percent/100,
         age.at.S1K = age) %>%
  select("mouseID", contains("S1K"), contains("total"))

GFP_data <- readxl::read_excel(file.path("datafiles", "RagGFP_ontogeny_pooled.xlsx"), sheet = 1)%>%
  select('mouseID', 'age', 'SP.8nai', 'LN.8nai') %>%
  mutate(total_counts = LN.8nai + SP.8nai,
         age.at.S1K = age) %>%
  select("mouseID", contains("S1K"), contains("total")) %>%
  left_join(GFPki_df, by = c('mouseID', 'age.at.S1K'))
  

ggplot() +
  geom_point(data=chimera_data, aes(x= age.at.S1K, y = total_counts), col=2)+
  geom_point(data=GFP_data, aes(x= age.at.S1K, y = total_counts), col=4) +
  scale_y_log10() + scale_x_log10()

ggplot() +
  geom_point(data=chimera_data, aes(x= age.at.S1K, y = total_kiprop*100), col=2)+
  geom_point(data=GFP_data, aes(x= age.at.S1K, y = total_kiprop*100), col=4) +
  scale_y_log10() + scale_x_log10()

## time points for solving the PDE
data_time_ont <- GFP_data$age.at.S1K                        # timepoints in observations 
data_time_chi <- chimera_data$age.at.S1K 
tb_chi <- chimera_data$age.at.BMT
ont_counts <- GFP_data$total_counts
ont_ki <- GFP_data$total_kiprop
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

for (n_shards in c(numObs)) {
  stan_rdump(c("numOnt",  "ont_ki", "ont_counts", 
               "numChi", "chi_counts",  "N_donor_fraction", "donor_ki", "host_ki",
               "ts_pred_ont", "ts_pred_chi1", "ts_pred_chi2", "ts_pred_chi3",
               "tb_pred1", "tb_pred2", "tb_pred3", "numPred", "n_shards", "dat_time", "dat_t0"),
             file = file.path('datafiles', paste0(Population, '_data_s', n_shards,".Rdump")))
} 

#### Model simulations
rstan::expose_stan_functions("stan_models/only_onto/MAP_asm_deltavar_cd4.stan")
#params <- c(0.3, 0.08, 0.005, 1e5, 3e5, 0.8)
params <- c(0.3, 0.04, 0.08, 0.005, 0.05, 0.01, 1e5, 0.8)
theta <- c(0)
x_i <- c(dat_t0[1])
x_r <- c(dat_time[1])
math_reduce(params, theta, 90, 40)

ts_pred_chi1 <- seq(66, 450, length.out = 300)
##
chi_vec <- sapply(ts_pred_chi1 - tb_pred1, Chi_spline)
total_counts <- N_pooled_time(ts_pred_chi1, tb_pred1, params)
donor_counts <- N_donor_time(ts_pred_chi1, tb_pred1, params)
Nfd <- donor_counts/(total_counts * chi_vec)

ggplot()+
  geom_line(aes(x = ts_pred_chi1, y=Nfd), size=2)


kseq <- seq(0, 1, 0.001)
ki_ini <- sapply(kseq, ki_dist_init)
ggplot()+
  geom_point(aes(x = kseq, y=ki_ini), size=2)

k_vec <- Vectorize(ki_dist_init)
integrate(k_vec, lower = 1/exp(1), upper = 1)

ki_theta <- sapply(kseq, ki_dist_theta, time=140)
ggplot()+
  geom_point(aes(x = kseq, y=ki_theta), size=2)

ageseq <- seq(0, 450)
age_dist <- sapply(ageseq, g_age, parms=params)
lambda <- sapply(ageseq, lambda_age, parms=params)
ggplot()+
  geom_point(aes(x = ageseq, y=age_dist), size=2)
ggplot()+
  geom_point(aes(x = ageseq, y=lambda), size=2)

theta_d <- sapply(ageseq, theta_donor, parms=params)
theta_h <- sapply(ageseq, theta_host, parms=params)
ggplot()+
  geom_point(aes(x = ageseq, y=theta_d), size=2)+
  geom_point(aes(x = ageseq, y=theta_h), col=2)


artf_dat <- data.frame()
artf_dat <- data.frame("x"= math_reduce(params, theta, data_time_chi[1], tb_chi[1]))
for (i in 2:length(data_time_chi)) {
  artf_dat[, i] = data.frame("x"=math_reduce(params, theta, data_time_chi[i], tb_chi[i]))
}
artf_df <- data.frame(t(artf_dat))

y3_mean <- c()
y4_mean <- c()
y5_mean <- c()
y6_mean <- c()

set.seed(1357)
for (i in 1:numChi){
  y3_mean[i] = artf_df[i, 1] #+ rnorm(1, 1e5, 1e5)
  y4_mean[i] = artf_df[i, 2] #+ rnorm(1, 0.01, 0.01)
  y5_mean[i] = artf_df[i, 3] #+ rnorm(1, 0.02, 0.02)
  y6_mean[i] = artf_df[i, 4] #+ rnorm(1, 0.04, 0.04)
}

ggplot()+
  geom_point(aes(x = dat_time[1:numChi], y = (y3_mean)))

ggplot() +
  geom_point(aes(x = dat_time[1:numChi], y = (y4_mean))) 

ggplot()+
  geom_point(aes(x = dat_time[1:numChi], y = (y5_mean))) +
  geom_point(aes(x = dat_time[1:numChi], y = (host_ki)), col=2)+ylim(0,1)

ggplot()+
  geom_point(aes(x = dat_time[1:numChi], y = (y6_mean))) +
  geom_point(aes(x = dat_time[1:numChi], y = (donor_ki)), col=2)+ylim(0,1)




pars <- c(out_table$mean[1:6])
chi_vec <- sapply(data_time_chi - tb_chi, Chi_spline)
init_conds <- c(872306.271*0.79226431, 872306.271 * (1-0.79226431), 0, 0)
yModel <- solve_ode_chi(data_time_chi, tb_chi, init_conds, pars)
stan_model_df <- data.frame("time" = data_time_chi,
                            "y_pred" = matrix(unlist(yModel), nrow = length(yModel), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         Nfd = (y_pred.1 + y_pred.2)/ (total_counts * chi_vec),
         host_ki = y_pred.3/(y_pred.1 + y_pred.4),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))

stan_model_df <- data.frame("time" = data_time_chi,
                            "y_pred" = matrix(unlist(yModel), nrow = length(yModel), byrow = TRUE))%>%
                              mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 +  y_pred.5 + y_pred.6 + y_pred.7 + y_pred.8,
                                     Nfd = (y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)/ (total_counts * chi_vec),
                                     host_ki = y_pred.3/(y_pred.3 + y_pred.4),
                                     donor_ki = y_pred.1/(y_pred.1 + y_pred.2))
logit_transf <- function(x) log(x/(1-x))
asinsqrt_transf <- function(x) asin(sqrt(x))

Resid_counts <- log(stan_model_df$total_counts) - log(chimera_data$total_counts)
Resid_Nfd <- logit_transf(stan_model_df$Nfd) - logit_transf(chimera_data$Nfd)
Resid_hostki <- asinsqrt_transf(stan_model_df$host_ki) - asinsqrt_transf(host_ki_df$ki_prop)
Resid_donorki <- asinsqrt_transf(stan_model_df$donor_ki) - asinsqrt_transf(donor_ki_df$ki_prop)

ggplot() +
  geom_point(aes(data_time_chi, y = Resid_counts), size=2, shape=1) +
  geom_hline(yintercept = 0.0, col ="darkred", size=1.0,  linetype=2)

ggplot() +
  geom_point(aes(data_time_chi, y = Resid_Nfd), size=2, shape=1) +
  geom_hline(yintercept = 0.0, col ="darkred", size=1.0,  linetype=2)

ggplot() +
  geom_point(aes(data_time_chi, y = Resid_hostki), size=2, shape=1) +
  geom_hline(yintercept = 0.0, col ="darkred", size=1.0,  linetype=2)

ggplot() +
  geom_point(aes(data_time_chi, y = Resid_donorki), size=2, shape=1) +
  geom_hline(yintercept = 0.0, col ="darkred", size=1.0,  linetype=2)

ypred_ont <- solve_ode_ont(data_time_ont, init_conds, pars)
stan_pred_ont <- data.frame("time" = data_time_ont,
                            "y_pred" = matrix(unlist(ypred_ont), nrow = length(ypred_ont), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         total_ki = (y_pred.1 + y_pred.3)/total_counts)

Resid_Ontcounts <- log(stan_pred_ont$total_counts) - log(ont_counts)
Resid_Ontki <- logit_transf(stan_pred_ont$total_ki) - logit_transf(ont_ki)

ggplot() +
  geom_point(aes(data_time_ont, y = Resid_Ontcounts), size=2, shape=1) +
  geom_hline(yintercept = 0.0, col ="darkred", size=1.0,  linetype=2)

ggplot() +
  geom_point(aes(data_time_ont, y = Resid_Ontki), size=2, shape=1) +
  geom_hline(yintercept = 0.0, col ="darkred", size=1.0,  linetype=2)


#ypred_ont <- solve_ode_ont(ts_pred_ont, init_conds, pars)
#stan_pred_ont <- data.frame("time" = ts_pred_ont,
#                            "y_pred" = matrix(unlist(ypred_ont), nrow = length(ypred_ont), byrow = TRUE))%>%
#  mutate(total_counts = y_pred.1 + y_pred.2,
#         total_ki = y_pred.1/(y_pred.1 + y_pred.2))
#
ypred1 <- solve_ode_chi(ts_pred_chi1, tb_pred1, init_conds, pars)
stan_pred_df1 <- data.frame("time" = ts_pred_chi1,
                            "y_pred" = matrix(unlist(ypred1), nrow = length(ypred1), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         fd = (y_pred.1 + y_pred.2)/ (total_counts * chi_vec),
         host_ki = y_pred.3/(y_pred.1 + y_pred.4),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))

ypred2 <- solve_ode_chi(ts_pred_chi2, tb_pred2, init_conds, pars)
stan_pred_df2 <- data.frame("time" = ts_pred_chi2,
                             "y_pred" = matrix(unlist(ypred2), nrow = length(ypred2), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         fd = (y_pred.1 + y_pred.2)/ (total_counts * chi_vec),
         host_ki = y_pred.3/(y_pred.1 + y_pred.4),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))

ypred3 <- solve_ode_chi(ts_pred_chi3, tb_pred3, init_conds, pars)
stan_pred_df3 <- data.frame("time" = ts_pred_chi3,
                            "y_pred" = matrix(unlist(ypred3), nrow = length(ypred3), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         fd = (y_pred.1 + y_pred.2)/ (total_counts * chi_vec),
         host_ki = y_pred.3/(y_pred.1 + y_pred.4),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))

###RTE
#ypred_ont <- solve_ode_ont(ts_pred_ont, init_conds, pars)
#stan_pred_ont <- data.frame("time" = ts_pred_ont,
#                            "y_pred" = matrix(unlist(ypred_ont), nrow = length(ypred_ont), byrow = TRUE))%>%
#  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
#         fd = (y_pred.1 + y_pred.2)/ (total_counts * chi_vec),
#         host_ki = y_pred.3/(y_pred.1 + y_pred.4),
#         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))
#
#
##init_conds <- c(0,0,0,0,874992*0.2, 874992 * 0.8, 0, 0)
#ypred1 <- solve_ode_chi(ts_pred_chi1, tb_pred1, init_conds, pars)
#stan_pred_df1 <- data.frame("time" = ts_pred_chi1,
#                            "y_pred" = matrix(unlist(ypred1), nrow = length(ypred1), byrow = TRUE))%>%
#  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 +  y_pred.5 + y_pred.6 + y_pred.7 + y_pred.8,
#         fd = (y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)/ (total_counts * chi_vec),
#         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
#         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))
#
#ypred2<- solve_ode_chi(ts_pred_chi2, tb_pred2, init_conds, pars)
#stan_pred_df2 <- data.frame("time" = ts_pred_chi2,
#                            "y_pred" = matrix(unlist(ypred2), nrow = length(ypred2), byrow = TRUE))%>%
#  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 +  y_pred.5 + y_pred.6 + y_pred.7 + y_pred.8,
#         fd = (y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)/ (total_counts * chi_vec),
#         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
#         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))
#
#ypred3<- solve_ode_chi(ts_pred_chi3, tb_pred3, init_conds, pars)
#stan_pred_df3 <- data.frame("time" = ts_pred_chi3,
#                            "y_pred" = matrix(unlist(ypred3), nrow = length(ypred3), byrow = TRUE))%>%
#  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 +  y_pred.5 + y_pred.6 + y_pred.7 + y_pred.8,
#         fd = (y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)/ (total_counts * chi_vec),
#         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
#         donor_ki = y_pred.1/(y_pred.1 + y_pred.2))
#
ggplot(stan_pred_df1)+
  geom_point(aes(x = time, y = (fd)))


ggplot()+
  geom_line(data = stan_pred_ont, aes(x = time, y = total_counts), col=2, size=1.5) +
  geom_line(data = stan_pred_df1, aes(x = time, y = total_counts), col=4, size=1.5)+
  geom_line(data = stan_pred_df2, aes(x = time, y = total_counts), col=3, size=1.5)+
  geom_point(data = stan_pred_df3, aes(x = time, y = total_counts), col=6) +
  scale_y_log10()


artf_dat <- data.frame()
artf_dat <- data.frame("x"= math_reduce(params, theta, data_time_chi[1], tb_chi[1]))
for (i in 2:length(data_time_chi)) {
  artf_dat[, i] = data.frame("x"=math_reduce(params, theta, data_time_chi[i], tb_chi[i]))
}
artf_df <- data.frame(t(artf_dat))

y1_mean <- c()
y2_mean <- c()
y3_mean <- c()
y4_mean <- c()
y5_mean <- c()
y6_mean <- c()

numOnt = 0
set.seed(1357)
for (i in 1:numOnt){
  y1_mean[i] = artf_df[i, 1] #+ rnorm(1, 1e5, 1e5)
  y2_mean[i] = artf_df[i, 2] #+ rnorm(1, 0.02, 0.02)
}

for (i in 1:numChi){
  y3_mean[i] = artf_df[numOnt + i, 1] #+ rnorm(1, 1e5, 1e5)
  y4_mean[i] = artf_df[numOnt + i, 2] #+ rnorm(1, 0.01, 0.01)
  y5_mean[i] = artf_df[numOnt + i, 3] #+ rnorm(1, 0.02, 0.02)
  y6_mean[i] = artf_df[numOnt + i, 4] #+ rnorm(1, 0.04, 0.04)
}


ggplot()+
  geom_point(aes(x = dat_time[1:numOnt], y = log(y1_mean)), size=3)+
  geom_point(aes(x = dat_time[(1+numOnt):(numChi+numOnt)], y = log(y3_mean)), col=2)

ggplot()+
  geom_point(aes(x = dat_time[1:numOnt], y = (y2_mean* 100)))+
  scale_y_log10(limits=c(0.3, 100)) +
  scale_x_log10()

ggplot()+
  geom_point(aes(x = dat_time[(1+numOnt):(numChi+numOnt)], y = log(y3_mean)))

ggplot()+
  geom_point(aes(x = dat_time[(1+numOnt):(numChi+numOnt)], y = (y4_mean)))

ggplot()+
  geom_point(aes(x = dat_time[(1+numOnt):(numChi+numOnt)], y = (y5_mean)))+
  ylim(0,1)

ggplot()+
  geom_point(aes(x = dat_time[(1+numOnt):(numChi+numOnt)], y = (y6_mean)))+
  ylim(0,1)

