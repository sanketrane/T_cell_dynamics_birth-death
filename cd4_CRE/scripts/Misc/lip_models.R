
library(readxl)
library(tidyverse)
library(deSolve)  

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
### density dependent proliferation (LIP) model
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
         group_id = ifelse(age_group == 1, 1,
                           ifelse(age_group == 7, 2,
                                  ifelse(age_group == 28, 3,
                                         ifelse(age_group == 56, 4, 5)))),
         time_track = age_anim - t0_group) %>%
  rename(mouse_id = id) 


### rate of loss -- changes with population density
lambda_dens <- function(Time, parms){
  lambda = parms[1]
  count_bar = parms[3]
  
  ## population density of the whole naive T cell compartment
  ## calculated using the spline
  pop_dens <- counts_func(Time, 7.06867, 5.54431, 0.15014)
  
  ### Net loss rate lambda = delta - rho
  value = lambda / (1 + (pop_dens/10^count_bar)^3)
  return(value)
}

pars_vec <- c(delta_log = -2, rho_log = -6, count_bar = 6)

#### fraction of ki67+ cells
ode_solver <- function(parms_vec, cohort_id, N0){
  filtered_df <- naive_df %>% 
    filter(group_id == cohort_id)
  
  ## time points for solving the PDE
  data_time <- filtered_df$time_track                        # time-points in observations 
  solve_time <- data_time %>% unique()          #unique time points to solve odes  
  time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
  
  ## ode function
  ode_func <-  function(Time, state, parms){
    with(as.list(c(state, parms)),{
      
      #labelled cells
      dY1 = -lambda_dens(Time, parms) * Y1
      
      #return the rate of change
      list(c(dY1))
      
    })  # end with(as.list ...
  }
  
  #initial conditions
  state <- c(10^N0)
  names(state) <- c("Y1")
  
  sol_ode <- ode(y=state, times=solve_time, func = ode_func, parms = parms_vec)
  
  return(sol_ode)
}


ode_solver(pars_vec, 6)

cohort1_preds <- 
naive_pred_df <- c(ode_solver(pars_vec, 1)[, 2],
                   ode_solver(pars_vec, 2)[, 2],
                   ode_solver(pars_vec, 3)[, 2],
                   ode_solver(pars_vec, 4)[, 2],
                   ode_solver(pars_vec, 5)[, 2])

naive_sim_df <- data.frame()
names(naive_sim_df) <- c("age_anim", "naive_labelled")
naive_pred <- naive_sim_df[time_index, 2]



optim_func <- function(parms, boot_data){
  
  ## predictions from the model
  fit_sim_df <- data.frame(ode_solver(solve_time, parms))
  names(fit_sim_df) <- c("age_anim", "naive_labelled")
  fit_pred <- fit_sim_df[time_index, 2]
  
  ssr = sum((log(boot_data$naive_labelled) - log(fit_pred))^2)
  
  k = length(parms)
  n1 = nrow(boot_data)
  
  logL = - (n1/2) * log(ssr)
  
  return(-logL)
}

counts_fit <- optim(par = pars_vec, fn = optim_func,
                     boot_data = naive_df, hessian = T, control = list(trace =T, maxit = 1000))

counts_fit$par


naive_sim_df <- data.frame(ode_solver(solve_time, counts_fit$par))
names(naive_sim_df) <- c("age_anim", "naive_labelled")
naive_pred <- naive_sim_df[time_index, 2]

ggplot()+
  geom_point(data=naive_df,  aes(x = age_anim, y = naive_labelled,
                                 col = as.factor(mouse_id)), alpha = 0.8, size =2) +
  geom_line(data=naive_sim_df, aes(x = age_anim, y = naive_labelled), size =1.2) +
  scale_y_log10(limits = c(5e3, 1e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 130), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  +
  myTheme + guides(color = F)



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
                                         ifelse(age_group == 56, N0_list$N0_g4, N0_list$N0_g5)))),
         group_id = ifelse(age_group == 1, 1,
                           ifelse(age_group == 7, 2,
                                  ifelse(age_group == 28, 3,
                                         ifelse(age_group == 56, 4, 5))))) %>%
  rename(mouse_id = id)


library(nlme)
datasetG <- groupedData(naive_labelled ~ 1 | mouse_id, naive_df)
nlin.mix <- nlme(fixed = log(naive_labelled) ~ log(n_total_vec(Time=age_anim, dstmp=age_group, t0=t0_group,
                                              )),
                 data = datasetG, 
                 random =  ~ sig_N0 | mouse_ID,
                 start = c("N0_mean" = 6, "l0"=0.05, "r_l"=0.01))

summary(nlin.mix)
fm2 <- update(nlin.mix, random = pdDiag(l0 + r_l ~ 1))
summary(fm2)

parms <- nlin.mix$coefficients$fixed
parms <- fm2$coefficients$fixed

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
  geom_line(data= filter(naive_df,  age_group ==7),  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(id)),size =2) +
  #scale_color_manual(values = c(1, 2, 3, 4, 5))+
  #scale_color_discrete(name = "Host age timestamp", label = c('d1', 'd7', 'd28', 'd56', 'd175'))+
  #geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled), size =1.2) +
  scale_y_log10(limits = c(5e3, 1e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 120), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  facet_wrap(~dstmp_group, scales = 'free_x') +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", colour = NA))

ggsave(filename = 'asm_a5dist_flat_p1.pdf', width = 8, height = 5, units = "in")

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

ggsave(filename = 'asm_a5dist_flat_p2.pdf', width = 6, height = 4.5, units = "in")








