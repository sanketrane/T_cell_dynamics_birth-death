setwd("~/Desktop/GIt_repos/magnum_opus")

library(readxl)
library(tidyverse)

theme_set(theme_bw())
myTheme <- theme(text = element_text(size = 11), axis.text = element_text(size = 11), axis.title =  element_text(size = 10, face = "bold"),
                 plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())


### importing data files
## RFP background correction
rfp_bg <- read_excel("data/timestamp/RFP_bg.xlsx", sheet = 1) %>%
  gather(key = age_cell, value = rfp_bg) %>% na.omit() %>%
  mutate(host_age = ifelse(age_cell == "d14", 14,
                           ifelse(age_cell == "d28", 28,
                                  ifelse(age_cell == "d56", 56,
                                         ifelse(age_cell == "d84", 84,245)))))

## function to give % of rfp background
rfp_func <- function(t, p, b0, b){
  
  p * b0 * exp(b*t)/(p + b0*(exp(b*t) - 1))
}

## fitting to get params
rfp_fit <- nls(rfp_bg ~ rfp_func(host_age, p, b0, b), data = rfp_bg,
               start = list(p = 0.2, b0 = 0.01, b = 0.01))

# results
summary(rfp_fit)


## Total counts of 
total_counts <- read_excel("data/timestamp/total_counts.xlsx", sheet = 1) %>%
  gather(key = age_cell, value = counts) %>% na.omit() %>%
  mutate(host_age = ifelse(age_cell == "d14", 14,
                           ifelse(age_cell == "d28", 28,
                                  ifelse(age_cell == "d56", 56,
                                         ifelse(age_cell == "d84", 84,245)))))
  
## function to give total counts of cd8 cells
counts_func <- function(t, k, n0, g){
  
    10^(k) * 10^(n0) * exp(g*t)/(10^(k) + 10^(n0) * (exp(g*t) - 1))
}


## fitting to get params
counts_fit <- nls(counts ~ counts_func(host_age, k, n0, g), data = total_counts,
               start = list(k = 7, n0 = 6, g = 0.1), control = list(maxiter = 500))

# results
summary(counts_fit)


## plots
ggplot()+
  geom_point(data = rfp_bg, aes(x = host_age, y = rfp_bg), size =2) + 
  geom_line(aes(x= ts_p, y = rfp_func(ts_p, 0.19455, 0.01647, 0.09148)), size =1.25) #+
  #scale_y_log10() #+ scale_x_log10(breaks = c(20, 60, 200, 600))

ggplot()+
  geom_point(data = total_counts, aes(x = host_age, y = counts), size =2) +
  geom_line(aes(x= ts_p, y = counts_func(ts_p, 7.06867, 5.54431, 0.15014)), size =1.25) +
  scale_y_log10(limits = c(5e5, 5e7)) + scale_x_log10(breaks = c(20, 60, 200, 600))



################################################################################
### Time stamp
tstmp <- read_excel("data/timestamp/data_03.xlsx", sheet = "data") 

tstmp_df <- tstmp %>%
  mutate(labelled_counts = (((fraction - rfp_func(age_anim, 0.19455, 0.01647, 0.09148))/100) 
                            * counts_func(age_anim, 7.06867, 5.54431, 0.15014)),
         naive_labelled = labelled_counts * np_fract/100)


### plots for whole pop
## labeled fraction
ggplot()+
  geom_point(data = tstmp_df,  aes(x = age_anim, y = fraction, col = as.factor(age_group)), alpha = 0.7, size =2) +
  scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100))+
  scale_color_discrete(name = "Age at timestamp (d)", guide = guide_legend(nrow = 3)) +
  labs(x = "Mouse age", title = "Counts of timestamped CD8 cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.88, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())

## labeled counts (corrected)
ggplot()+
  geom_point(data = tstmp_df,  aes(x = age_anim, y = labelled_counts, col = as.factor(age_group)), size =2) +
  scale_y_log10(limits = c(1e4, 1e7), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  scale_color_discrete(name = "Age at timestamp (d)", guide = guide_legend(nrow = 3)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped CD8 T cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.88, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())


## plots for naive pop
## labeled counts (corrected)
ggplot()+
  geom_point(data = tstmp_df,  aes(x = age_anim, y = naive_labelled,
                                   col = as.factor(age_group)), alpha = 0.8, size =2.2) +
  scale_y_log10(limits = c(5e3, 2.5e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  scale_color_discrete(name = "Age at timestamp (d)", guide = guide_legend(nrow = 3)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.85, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())


################################################################################
### Age structured model

# initial age distribution at host age t0 == t_init
g_age <- function(age, dstmp, t0, N0){
  if(age < t0 - dstmp - 5){
   value = 0.0
  } else if(age <= t0 - dstmp) {
   value = (10^N0)/5
  } else{
    value = 0.0
  }
}
g_age_vec <- Vectorize(g_age)


## rate of loss -- changes with cell age
lambda_age <- function(age, l0, r_l){
  value = l0 * exp(- r_l * age)
  return(value)
}

## total cell age distribution atr time == Time
G_age_time <- function(age, Time, dstmp, t0, N0, l0, r_l){
  g_age_vec(age - Time + t0, dstmp, t0, N0) *
           exp(- integrate(lambda_age, lower=age - Time + t0, upper=age, l0=l0, r_l=r_l)$value)
}

## vectoring function
G_age_vec <- Vectorize(G_age_time)

## total counts of labelled cells at time == Time
n_total_time <- function(Time, dstmp, t0, N0, l0, r_l){
    integrate(G_age_vec, lower = Time - dstmp - 5, upper = Time - dstmp, Time=Time,
            dstmp=dstmp, t0=t0, N0=N0, l0=l0, r_l=r_l)$value
}

## vectoring function
n_total_vec <- Vectorize(n_total_time)

##  function used by the original paper -- exponential decay 
labelled_func <- function(cell_age, host_age, parms){
    N0 = parms[1]
  return(10^(N0) * exp(- lambda_age(cell_age, parms) * host_age))
}


### N0 for each cohort
N0_list <- list(
  N0_g1 = 5.24,
  N0_g2 = 5.88,
  N0_g3 = 5.81,
  N0_g4 = 5.27,
  N0_g5 = 5.29
)

### data for fitting
naive_df <- tstmp_df %>%
  select('id', 'age_anim', 'age_group', 'naive_labelled') %>%
  mutate(t0_group = ifelse(age_group == 1, 14,
                           ifelse(age_group == 7, 28,
                                  ifelse(age_group == 28, 42,
                                         ifelse(age_group == 56, 70, 189)))),
         dstmp_group = ifelse(age_group == 1, "Group1",
                              ifelse(age_group == 7, "Group2",
                                     ifelse(age_group == 28, "Group3",
                                            ifelse(age_group == 56, "Group4", "Group5")))),
         N0_group = ifelse(age_group == 1, N0_list$N0_g1,
                           ifelse(age_group == 7, N0_list$N0_g2,
                                  ifelse(age_group == 28, N0_list$N0_g3,
                                         ifelse(age_group == 56, N0_list$N0_g4, N0_list$N0_g5))))) 



## predictions for individual age cohorts
tseq1 <- seq(14, 112, length.out = 50)
sim1 <- n_total_vec(Time=tseq1, dstmp=1, t0=14, N0=N0_list$N0_g1, l0=0.07, r_l=0.018)
tseq2 <- seq(28, 112, length.out = 50)
sim2 <- n_total_vec(Time=tseq2, dstmp=7, t0=28, N0=N0_list$N0_g2, l0=0.07, r_l=0.018)
tseq3 <- seq(42, 98, length.out = 50)
sim3 <- n_total_vec(Time=tseq3, dstmp=28, t0=42, N0=N0_list$N0_g3, l0=0.07, r_l=0.018)
tseq4 <- seq(70, 154, length.out = 50)
sim4 <- n_total_vec(Time=tseq4, dstmp=56, t0=70, N0=N0_list$N0_g4, l0=0.07, r_l=0.018)
tseq5 <- seq(189, 273, length.out = 50)
sim5 <- n_total_vec(Time=tseq5, dstmp=175, t0=189, N0=N0_list$N0_g5, l0=0.07, r_l=0.018)


## labels for individual age cohorts
age_cohorts <- c(rep("Group1", 50), rep("Group2", 50), rep("Group3", 50), 
                 rep("Group4", 50), rep("Group5", 50))


## predictions data frame
naive_sim_df <- data.frame("age_anim" = c(tseq1, tseq2, tseq3, tseq4, tseq5),
                           "dstmp_group" = age_cohorts,
                           "naive_labelled" = c(sim1, sim2, sim3, sim4, sim5))



## plots for naive pop
## labeled counts (corrected) with predictions from ASM
ggplot()+
  geom_point(data = naive_df,  aes(x = age_anim, y = naive_labelled,
                                   col = dstmp_group),size =2) +
  scale_color_discrete(name = "Host age timestamp", label = c('d1', 'd7', 'd28', 'd56', 'd175'))+
  geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled), size =1.2) +
  scale_y_log10(limits = c(5e3, 2.5e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  facet_wrap(~dstmp_group, scales = 'free_x') +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", colour = NA))
  

ggplot()+
  geom_point(data = naive_df,  aes(x = age_anim, y = naive_labelled,
                                            col = dstmp_group), alpha = 0.8, size =2) +
  geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled, 
                                     col = dstmp_group), size =1.2) +
  scale_color_discrete(name = "Host age at timestamp", label = c('d1', 'd7', 'd28', 'd56', 'd175'),
                       guide = guide_legend(nrow = 2))+
  scale_y_log10(limits = c(5e3, 2.5e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.85, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())


### data for fitting
naive_fit_df <- naive_df %>%
  mutate(naive_fit = n_total_vec(Time=age_anim, dstmp=age_group, t0=t0_group, N0=N0_group, l0=0.07, r_l=0.018))

ggplot()+
  geom_point(data = naive_df,  aes(x = age_anim, y = naive_labelled,
                                   col = dstmp_group),size =2) +
  scale_color_discrete(name = "Host age timestamp", label = c('d1', 'd7', 'd28', 'd56', 'd175'))+
  geom_point(data = naive_fit_df, aes(x = age_anim, y = naive_fit), shape = 17, size =2.5) +
  scale_y_log10(limits = c(5e3, 2.5e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  facet_wrap(~dstmp_group, scales = 'free_x') +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", colour = NA))

ggplot()+
  geom_point(data = naive_df,  aes(x = age_anim, y = naive_labelled,
                                   col = dstmp_group), alpha = 0.8, size =2) +
  geom_point(data = naive_fit_df, aes(x = age_anim, y = naive_fit,
                                      col = dstmp_group), shape = 17, size =4) +
  scale_color_discrete(name = "Host age at timestamp", label = c('d1', 'd7', 'd28', 'd56', 'd175'),
                       guide = guide_legend(nrow = 2))+
  scale_y_log10(limits = c(5e3, 2.5e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.85, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())


#fitting the ASM to log(labeled naive counts)
LL_naive_labelled <- function(param, boot_data) { 
  l0  <- param[1]
  r_l  <- param[2]             #parameters to be estimated as part of a vector
  
  k  <- length(param)          #number of unknown parameters 
  n1 <- nrow(boot_data)        #number of observations in the data-set
  
  fit_df <- boot_data %>%
    mutate(naive_fit = n_total_vec(Time=age_anim, dstmp=age_group, t0=t0_group, N0=N0_group, l0=l0, r_l=r_l))
  
  ssr1 <- sum((log(boot_data$naive_labelled) - log(fit_df$naive_fit))^2)   #SSR for the data-set
  
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multiplication of residual and transpose of residuals
  logl <- -(n1/2) * log(ssr1) 
  
  aiccd8 <<- -2*logl + 2*k
  
  return(-logl)     #since optim minimizes the function by default, ML
} 

fit_naive_labelled <- optim(par=c(0.05,0.01), fn=LL_naive_labelled, boot_data=naive_df,
                            control = list(trace = TRUE, maxit = 1000))
fit_naive_labelled
aiccd8
parms <- fit_naive_labelled$par

## predictions for individual age cohorts
tseq1 <- seq(10, 120, length.out = 50)
sim1 <- n_total_vec(Time=tseq1, dstmp=1, t0=14, N0=N0_list$N0_g1, l0=parms[1], r_l=parms[2])
tseq2 <- seq(20, 120, length.out = 50)
sim2 <- n_total_vec(Time=tseq2, dstmp=7, t0=28, N0=N0_list$N0_g2, l0=parms[1], r_l=parms[2])
tseq3 <- seq(35, 105, length.out = 50)
sim3 <- n_total_vec(Time=tseq3, dstmp=28, t0=42, N0=N0_list$N0_g3, l0=parms[1], r_l=parms[2])
tseq4 <- seq(60, 160, length.out = 50)
sim4 <- n_total_vec(Time=tseq4, dstmp=56, t0=70, N0=N0_list$N0_g4, l0=parms[1], r_l=parms[2])
tseq5 <- seq(180, 280, length.out = 50)
sim5 <- n_total_vec(Time=tseq5, dstmp=175, t0=189, N0=N0_list$N0_g5, l0=parms[1], r_l=parms[2])


## labels for individual age cohorts
age_cohorts <- c(rep("Group1", 50), rep("Group2", 50), rep("Group3", 50), 
                 rep("Group4", 50), rep("Group5", 50))


## predictions data frame
naive_sim_df <- data.frame("age_anim" = c(tseq1, tseq2, tseq3, tseq4, tseq5),
                           "dstmp_group" = age_cohorts,
                           "naive_labelled" = c(sim1, sim2, sim3, sim4, sim5))
## plots for naive pop
## labeled counts (corrected) with predictions from ASM
ggplot()+
  geom_point(data=naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = dstmp_group),size =2) +
  scale_color_discrete(name = "Host age timestamp", label = c('d1', 'd7', 'd28', 'd56', 'd175'))+
  geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled), size =1.2) +
  scale_y_log10(limits = c(5e3, 2.5e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  facet_wrap(~dstmp_group, scales = 'free_x') +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", colour = NA))


ggplot()+
  geom_point(data=naive_df,  aes(x = age_anim, y = naive_labelled,
                                   col = dstmp_group), alpha = 0.8, size =2) +
  geom_line(data=naive_sim_df, aes(x = age_anim, y = naive_labelled, 
                                     col = dstmp_group), size =1.2) +
  scale_color_discrete(name = "Host age at timestamp", label = c('d1', 'd7', 'd28', 'd56', 'd175'),
                       guide = guide_legend(nrow = 2))+
  scale_y_log10(limits = c(5e3, 2.5e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.85, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())










