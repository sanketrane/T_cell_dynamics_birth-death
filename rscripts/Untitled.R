rm(list = ls()); gc():
  
library(deSolve)
library(tidyverse)
library(rstan)


rstan::expose_stan_functions("stan_models/only_chimera/MAP_rtemld_cd4.stan")

Ode_solver <- function(time_vec, init_cond, parms_vec){
  
  ode_func <-  function(t, state, parms){
    with(as.list(c(state, parms)),{
      
      delta_rte = parms[1]
      rho_rte = parms[2]
      delta_nai = parms[3]
      rho_nai = parms[4]
      mu = parms[5]
      
      beta = 1/3.5;
      
      #RTE Ki67hi
      dY1 = rho_rte * (1 - mu) * (2 * Y2 + Y1) - (beta + delta_rte + (rho_rte * mu)) * Y1
      #RTE Ki67lo
      dY2 = beta * Y1 - (rho_rte + delta_rte) * Y2
      # mN Ki67hi
      dY3 = 2 * rho_rte * mu * (Y2 + Y1) + rho_nai * (2 * Y4 + Y3) - (beta + delta_nai) * Y3;
      # mN Ki67lo
      dY4 = beta * Y3 - (rho_nai + delta_nai) * Y4;  
      
      #return the rate of change
      list(c(dY1, dY2, dY3, dY4))
      
    })  # end with(as.list ...
  }
  
  # time points for which conc is reported
  # include the points where data is available
  times = time_vec
  sol_ode <- ode(y= init_cond, times=times, func = ode_func, parms = parms_vec)
  
  sol_df <-  data.frame(sol_ode) %>%
    mutate(total_counts = Y1 + Y2 + Y3 + Y4) %>% 
    select("time", "total_counts")
  
  return(sol_df$total_counts)
}


## prediction points
tseq1 <- seq(14, 300, length.out = 200)
tseq2 <- seq(28, 300, length.out = 200)
tseq3 <- seq(42, 300, length.out = 200)
tseq4 <- seq(70, 300, length.out = 200)
tseq5 <- seq(189, 300, length.out = 200)

td_df <- data.frame(tseq1, tseq2, tseq3, tseq4, tseq5)

#initial conditions
par_est <- c(0.698799621, 0.004285617, 0.033163141, 0.000479521, 0.001063894, 0.087248909)
data_time <- c(1, 14, 28, 42, 70, 189)
init_conds <- c(897398.161*0.802040798, 897398.161 * (1-0.802040798), 0, 0)
yModel <- solve_ode_ont(data_time, init_conds, par_est)
stan_pred_df <- data.frame("time" = data_time,
                            "y_pred" = matrix(unlist(yModel), nrow = length(yModel), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)



par_vec <-  c(0.033163141, 0.001063894, 0.004285617, 0.000479521, 0.087248909)
counts_df <- data.frame()
for (i in 1:5) {
  init_transfer <- c("Y1" = stan_pred_df$y_pred.1[i+1], 
                 "Y2" = stan_pred_df$y_pred.2[i+1], 
                 "Y3" = stan_pred_df$y_pred.3[i+1],  
                 "Y4" = stan_pred_df$y_pred.4[i+1])
  counts_df[(1 + (i - 1) * 200): (i * 200), 'Time_since_treatment'] <- td_df[1:200, i] #- data_time[i+1]
  counts_df[(1 + (i - 1) * 200): (i * 200), 'Host_age'] <- data_time[i+1]
  counts_df[(1 + (i - 1) * 200): (i * 200), 'Cell_counts'] <- Ode_solver(td_df[1:200, i] - data_time[i+1], init_transfer, par_vec) 
  
}


ggplot(counts_df)+
  geom_point(aes(x = Time_since_treatment, y = Cell_counts, col = as.factor(Host_age)), size = 1.2) +
  scale_y_log10(limits = c(5e3, 3e7)) +
 labs(x = 'Host age', y = NULL, title = paste0("Cell counts of RFP cohorts")) +
  #scale_x_log10(limits = c(1, 450), breaks = c(10, 30, 100, 300, 450)) +
  theme(legend.position = c(0.15, 0.1))














