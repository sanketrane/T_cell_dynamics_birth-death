library(tidyverse)
library(deSolve)

setwd("~/Desktop/GIt_repos/magnum_opus")

### data import

dat_cd4 <- read.csv("data/naive_cd8.csv") %>% arrange(time)

# unique time points in data for odes solver
unique_times_df <- dat_cd4 %>% distinct(time, .keep_all = TRUE)

data_time <- dat_cd4$time                                          # data and solver time 
solve_time <- unique_times_df$time                                       # unique time points to solve ode
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time

#fixed parameters
Beta = 1/3.5
Beta_prime = 1/2.5
ts_pred <- seq(3, 300, 0.5)

### transformation functions
logit_trans <- function(x){
  
  log(x/(1-x))
}

expit_trans <- function(x){
  
  exp(x)/(1 + exp(x))
}

#### total counts

# the function for total thymic output
theta_func <- function(Time, psi_logit) {
  
  psi = expit_trans(psi_logit)
  
  basl  = 2.5e5
  theta = 2e3
  n     = 2
  X     = 34
  q     = 3.5
  
  theta = basl + (theta * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q))))
  
  return(psi * theta)
}

eps_func <- function(Time){
  
  eps = exp(-0.03337899 * (Time + 2.92554110)) + 0.13732103
  
  return(eps)
}

#ggplot() + geom_point(aes(ts_pred, theta_func(ts_pred, -0.85)))

### rate of turnover varying with cell age
delta_age <- function(age, delta0_log, r_log){
  
  ### transforming the parameters
  delta0 = exp(delta0_log)
  r = exp(r_log)
  
  delta0 * exp(-r * age)
}

#ggplot() + geom_point(aes(ts_pred, delta_age(ts_pred, -3, -5)))

lambda_age <- function(age, delta0_log, r_log, rho_log){
  
  ### transforming the parameters
  rho = exp(rho_log)
  
  delta_age(age, delta0_log, r_log) -  rho
}

#ggplot() + geom_point(aes(ts_pred, lambda_age(ts_pred, -3, -5, -5)))

# age distribution of cells
G_noabmt <- function(a, t, psi_logit,  delta0_log, r_log, rho_log) {  
  
  theta_func(t-a, psi_logit)  * exp(- integrate(lambda_age, lower = 0, upper = a, delta0_log = delta0_log, r_log = r_log, rho_log = rho_log)$value)
}  

G_vec <-  Vectorize(G_noabmt)

#ggplot() + geom_point(aes(ts_pred, G_vec(ts_pred, 600, 0.3, -3, -5, -5)))


N_G_noabmt <- function(t, psi_logit, delta0_log, r_log, rho_log){
  
  integrate(G_vec, lower = 0, upper = t, t = t, psi_logit = psi_logit, delta0_log = delta0_log, r_log = r_log, rho_log = rho_log)$value
}

N_G.v <- Vectorize(N_G_noabmt)

#ggplot() + geom_point(aes(ts_pred, N_G.v(ts_pred, -0.85, -3, -5, -5)))

counts_func <- function(ts, param){
  
  psi_logit  <- param[1]
  delta0_log <- param[2]                     #parametrs to be estimated as part of a vector
  r_log      <- param[3]        
  rho_log    <- param[4]        
  
  
  counts = N_G.v(ts, psi_logit, delta0_log, r_log, rho_log)
  
  return(counts)
}

##### Ki67 proportions

numer <-  function(age, Time, psi_logit, delta0_log, r_log, rho_log){
  
  delta_age(age, delta0_log, r_log) * G_vec(age, Time, psi_logit, delta0_log, r_log, rho_log)
}


## Delta(t) is the integral of delta(age) from t0, t
Delta <- function(Time, psi_logit, delta0_log, r_log, rho_log){
  
  ifelse(Time >0,
  (integrate(numer, lower = 0, upper = Time, Time = Time, psi_logit = psi_logit, delta0_log = delta0_log, r_log = r_log, rho_log = rho_log)$value)/
    N_G.v(Time, psi_logit, delta0_log, r_log, rho_log),
  exp(delta0_log))
}

Delta_vec <- Vectorize(Delta)

#ggplot() + geom_point(aes(ts_pred, Delta_vec(ts_pred, -0.85, -3, -5, -5)))


#### fraction of ki67+ cells
ki_func <- function(ts, param){
  
  ## ode function
   ode_func <-  function(Time, state, parms){
     with(as.list(c(state, parms)),{
       
       #Ki67 hi VRTE
       dY1 <- theta_func(Time, psi_logit) * eps_func(Time) - (Beta_prime + exp(rho_log) + Delta_vec(Time, psi_logit, delta0_log, r_log, rho_log)) * Y1
       
       #Ki67 hi mN
       dY2 <- exp(rho_log) * (2 * Y1 + Y2 + 2 * Y3) - (Beta + Delta_vec(Time, psi_logit, delta0_log, r_log, rho_log)) * Y2
       
       #Ki67 lo mN
       dY3 <- theta_func(Time, psi_logit) * (1 - eps_func(Time)) + Beta_prime * Y1 + Beta * Y2 - (exp(rho_log) + Delta_vec(Time, psi_logit, delta0_log, r_log, rho_log)) * Y3
       
       #return the rate of change
       list(c(dY1, dY2, dY3))
       
     })  # end with(as.list ...
   }
   
   #initial conditions
   #state <- c(X = N0 * (1 - eps), Y = N0 * eps)
   state <-  c(Y1 = 0, Y2 = 0, Y3 = 0)
   
   # time points for which conc is reported
   # include the points where data is available
   times = c(0, ts)
   
   sol_ode <- ode(y= state, times=times, func = ode_func, parms = param)
   
   sol_df <-  data.frame(sol_ode) %>%
     mutate(ki_prop = (Y1 + Y2)/ (Y1 + Y2 + Y3)) %>% na.omit
   
   return(sol_df$ki_prop)
}

#ggplot() + geom_point(aes(ts_pred, ki_func(ts_pred, pars)))


#fitting log N_total & logit N_fd to the NCD4 data
ADT_NCD4 <- function(param, boot_data) { 
  
  obs_data.df <- boot_data%>%
    mutate(log_counts = log(counts),
           logit_ki = logit_trans(ki67)) 
  
  pred_counts <- log(counts_func(solve_time, param))
  
  pred_ki <- logit_trans(ki_func(solve_time, param))
  
  ssqres1 <- sum((pred_counts[time_index] - obs_data.df$log_counts)^2)
  ssqres2 <- sum((pred_ki[time_index] - obs_data.df$logit_ki)^2)
  
  
  k  <- length(param)                #number of unknown parameters 
  n1 <- nrow(obs_data.df)              #number of observations in dataset1
  
  #cost function
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multipltication of residual and transpose of residuals
  logl <-  - (n1/2)* log(ssqres1)  - (n1/2)* log(ssqres2)  
  
  #stats for model validation
  AICcd4 <<-  2*k - 2*logl
  
  return(-logl)     #since optim minimizes the function by default, ML
} 

pars <- list(psi_logit = -0.5, delta0_log = -3, r_log = -5, rho_log = -3)

ADT_NCD4_fit <- optim(par = pars, fn = ADT_NCD4, boot_data = dat_cd4, hessian = TRUE, control = list(trace = T))
expit_trans(ADT_NCD4_fit$par[1])
1/exp(ADT_NCD4_fit$par[2])
log(2)/exp(ADT_NCD4_fit$par[3])
1/exp(ADT_NCD4_fit$par[4])
AICcd4

par_est <- c(ADT_NCD4_fit$par[1],  ADT_NCD4_fit$par[2], ADT_NCD4_fit$par[3], ADT_NCD4_fit$par[4])

pred_df <- data.frame("Host_age" = ts_pred,
                      "counts" = counts_func(ts_pred, par_est),
                      "ki_prop" = ki_func(ts_pred, par_est))


ggplot() +
  geom_point(data = dat_cd4, aes(time, y= counts), size = 2) +
  labs(title = "Naive CD4 cell counts", y=NULL, x="mouse age (days)") + 
  geom_line(data = pred_df, aes(Host_age, y= counts), color="#2548b4", size = 1.5)+
  scale_y_continuous(limits=c(1e5, 5e7),  trans="log10", breaks=c(1e7, 1e6, 1e8)) +  #xlim(0, 300) +
  scale_x_continuous(limits=c(3, 320),  trans="log10", breaks=c(3, 10, 30, 100, 300)) +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title =  element_text(size = 12, face = "bold"),
                     panel.background = element_rect(colour = "black"),
                     plot.title = element_text(size=12, face = "bold",  hjust = 0.5),
                     legend.text = element_text(size=12, hjust = 0.5),
                     legend.title = element_blank())
ggplot() +
  geom_point(data = dat_cd4, aes(time, y= ki67), size = 2) +
  labs(title="Ki67high proportions", y=NULL, x="mouse age (days)") +
  geom_line(data = pred_df, aes(Host_age, y= ki_prop), color="#2548b4", size = 1.5)+
  scale_x_continuous(limits=c(3, 320),  trans="log10", breaks=c(3, 10, 30, 100, 300)) +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title =  element_text(size = 12, face = "bold"),
                     panel.background = element_rect(colour = "black"),
                     plot.title = element_text(size=12, face = "bold",  hjust = 0.5),
                     legend.text = element_text(size=12, hjust = 0.5),
                     legend.title = element_blank())

ggplot()+
  geom_line(aes(x = ts_pred, y = delta_age(ts_pred, ADT_NCD4_fit$par[2], ADT_NCD4_fit$par[3])))



## open graphics device 
## saving  plots for quality control 
pdf(file = file.path("~/Dropbox", paste("ASM_onto","Plots%03d.pdf", sep = "")),
    width = 6, height = 4, onefile = F, useDingbats = FALSE)

#pairs(fit, pars = parametersToPlot)

p1 <- ggplot() +
  geom_point(data = dat_cd4, aes(time, y= counts), size = 2) +
  labs(title = "Naive CD4 cell counts", y=NULL, x="mouse age (days)") + 
  geom_line(data = pred_df, aes(Host_age, y= counts), color="#2548b4", size = 1.5)+
  scale_y_continuous(limits=c(1e5, 5e7),  trans="log10", breaks=c(1e7, 1e6, 1e8)) +  #xlim(0, 300) +
  scale_x_continuous(limits=c(3, 320),  trans="log10", breaks=c(3, 10, 30, 100, 300)) +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title =  element_text(size = 12, face = "bold"),
                     panel.background = element_rect(colour = "black"),
                     plot.title = element_text(size=12, face = "bold",  hjust = 0.5),
                     legend.text = element_text(size=12, hjust = 0.5),
                     legend.title = element_blank())


p2 <- ggplot() +
  geom_point(data = dat_cd4, aes(time, y= ki67), size = 2) +
  labs(title="Ki67high proportions", y=NULL, x="mouse age (days)") +
  geom_line(data = pred_df, aes(Host_age, y= ki_prop), color="#2548b4", size = 1.5)+
  scale_x_continuous(limits=c(3, 320),  trans="log10", breaks=c(3, 10, 30, 100, 300)) +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title =  element_text(size = 12, face = "bold"),
                     panel.background = element_rect(colour = "black"),
                     plot.title = element_text(size=12, face = "bold",  hjust = 0.5),
                     legend.text = element_text(size=12, hjust = 0.5),
                     legend.title = element_blank())

p3 <- ggplot()+
  geom_line(aes(x = ts_pred, y = Delta(ts_pred, ADT_NCD4_fit$par[1],  ADT_NCD4_fit$par[2], ADT_NCD4_fit$par[3], ADT_NCD4_fit$par[4])))


p1
p2
p3

dev.off()























