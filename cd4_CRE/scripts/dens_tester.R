rm(list = ls()); gc()
library(rstan)
library(tidyverse)

## importing ASM stan model
ModelName <- 'DDM_N0mouse'
stan_file <- paste0(ModelName, '.stan')
expose_stan_functions(file.path('stan_models', stan_file))

parms_vec <- c(0.4, 3.2e6)

#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv")

solve_time <- naive_df$age_anim %>% unique()

N0_mouse <- c(2e5, 1.5e5, 3e5)

N_total(naive_df$age_anim[2], naive_df$t0_group[2], N0_mouse[1], parms_vec)
N_total_time(solve_time, 2e5, parms_vec)


solve_time <- seq(14, 120, length.out = 100)

naive_sim_df <- data.frame("age_anim" =solve_time)
for (i in 1:length(solve_time)){
  naive_sim_df[i, "naive_labelled"] = solve_ode(solve_time, 2e5, parms_vec)[[i]]
}

total_counts_vec <- sapply(solve_time, counts_func)
 
ggplot()+
  geom_point(data=naive_df,  aes(x = age_anim, y = naive_labelled,
                                 col = as.factor(mouse_id)), alpha = 0.8, size =2) +
  geom_line(data=Y1pred, aes(x = timeseries, y = naive_labelled), size =1.2) +
  scale_y_log10(limits = c(5e3, 1e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 130), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  +
  guides(color = F)



### parameters table
num_pars <- 71 

parameters_to_plot <- model_fit@model_pars[1:8]

ptable <- monitor(as.array(model_fit, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]

## stan solution
## prediction points
tseq1 <- seq(14, 120, length.out = 200)
tseq2 <- seq(28, 120, length.out = 200)
tseq3 <- seq(42, 110, length.out = 200)
tseq4 <- seq(70, 160, length.out = 200)
tseq5 <- seq(189, 280, length.out = 200)

sim1 <- sapply(tseq1,  N_total, t0 = 14, N0 = mean(out_table$mean[1:12]),
               parms = c(out_table$mean[67], out_table$mean[72]))
sim2 <- sapply(tseq2,  N_total, t0 = 28, N0 = mean(out_table$mean[13:29]),
               parms = c(out_table$mean[68], out_table$mean[72]))
sim3 <- sapply(tseq3,  N_total, t0 = 42, N0 = mean(out_table$mean[30:41]),
               parms = c(out_table$mean[69], out_table$mean[72]))
sim4 <- sapply(tseq4,  N_total, t0 = 70, N0 = mean(out_table$mean[42:53]),
               parms = c(out_table$mean[70], out_table$mean[72]))
sim5 <- sapply(tseq5,  N_total, t0 = 189, N0 = mean(out_table$mean[54:66]),
               parms = c(out_table$mean[71], out_table$mean[72]))


naive_sim_df <- rbind(Y1pred, Y2pred, Y3pred, Y4pred, Y5pred)
naive_df$group_id <- as.factor(naive_df$group_id)
levels(naive_df$group_id) <- seq(1, 5)#c('Cohort-1', "Cohort-2", 'Cohort-3', "Cohort-4",'Cohort-5')

## labels for individual age cohorts
age_cohorts <- c(rep(1, 200), rep(2, 200), rep(3, 200), 
                 rep(4, 200), rep(5, 200))
## predictions data frame
naive_sim_df <- data.frame("age_anim" = c(tseq1, tseq2, tseq3, tseq4, tseq5),
                           "group_id" = age_cohorts,
                           "naive_labelled" = c(sim1, sim2, sim3, sim4, sim5))
ggplot()+
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(group_id)), size =3) +
  geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled, 
                                     col = as.factor(group_id)), size =1.2) +
  scale_color_discrete(name = "Cohort", guide = guide_legend(nrow = 3)) +
  scale_y_log10(limits = c(5e3, 2e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  +# myTheme +
  #facet_wrap(~ as.factor(group_id), scales = 'free_x') + guides(col = F) +
  theme(legend.position = c(0.9, 0.85),
        legend.background = element_rect(fill = NA, colour = NA))


ggplot()+
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(group_id)), size =3) +
  geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled, 
                                     col = as.factor(group_id)), size =2) +
  scale_color_discrete(name = "Cohort", guide = guide_legend(nrow = 3)) +
  scale_y_log10(limits = c(5e3, 2e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  +# myTheme +
  facet_wrap(~ as.factor(group_id), scales = 'free_x') + guides(col = F) 








###### ------- R solution ------- ########
counts_funcR <- function(t, k, n0, g){
  
  10^(k) * 10^(n0) * exp(g*t)/(10^(k) + 10^(n0) * (exp(g*t) - 1))
}

lambda_densR <- function(Time, parms){
  lambda = parms[1]
  count_bar = parms[2]
  
  ## population density of the whole naive T cell compartment
  ## calculated using the spline
  pop_dens <- counts_funcR(Time, 7.06867, 5.54431, 0.15014)
  
  ### Net loss rate lambda = delta - rho
  value = lambda / (1 + (pop_dens/count_bar)^3)
  return(value)
}

lambda_densR(c(14, 28), parms_vec)

## ode function
ode_func <-  function(Time, state, parms){
  with(as.list(c(state, parms)),{
    
    #labelled cells
    dY1 = -lambda_densR(Time, parms) * Y1
    
    #return the rate of change
    list(c(dY1))
    
  })  # end with(as.list ...
}

#initial conditions
state <- c(2e5)
names(state) <- c("Y1")

sol_ode <- ode(y=state, times=c(14, 28), func = ode_func, parms = parms_vec)
