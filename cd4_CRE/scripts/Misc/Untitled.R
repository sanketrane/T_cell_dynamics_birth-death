library(tidyverse)
library(deSolve)


#Known parameters
v_cd4 = 0.002146831
v_cd8 = 0.003014382
ta = 49

# peripheral turnover rate 
lambda_hill <- function(a,l,r) {
  l/(1+ (a/r)^2)
}

#total thymic output
theta_func <- function(Time, theta0) {
  
  return(10^theta0 * exp(-v_cd4 * Time))
}

### rate of turnover varying with cell age
delta_age <- function(age, delta0, r){
  
  delta0 * exp(-r * age)
}

lambda_age <- function(age, delta0, r, mu, rho){
  
  delta_age(age, delta0, r) + mu -  rho
}

# age distribution of cells
G_noabmt <- function(a, t, theta0, delta0, r, mu, rho) {  
  
  theta_func(t-a, theta0)  * exp(- integrate(lambda_age, lower = 0, upper = t, delta0 = delta0, r = r, mu = mu, rho = rho)$value)
}  

G_vec <-  Vectorize(G_noabmt)

G_vec(100, 100, 6, 0.02, 0.001, 0.01, 0.01)

N_G_noabmt <- function(t, theta0, delta0, r, mu, rho){
  
  integrate(G_vec, lower = 0, upper = t, t = t, theta0 = theta0, delta0 = delta0, r = r, mu = mu, rho = rho)$value
}

N_G.v <- Vectorize(N_G_noabmt)

N_G.v(100, 6, 0.02, 0.001, 0.01, 0.01)

numerator_delta <-  function(a, t,theta0, delta0, r, mu, rho){
  
  (delta_age(a, delta0, r) * G_vec(a, t, theta0, delta0, r, mu, rho))
}

normalised_death <- function(a, t,theta0, delta0, r, mu, rho){
  
  (delta_age(a, delta0, r) * G_vec(a, t, theta0, delta0, r, mu, rho))/
    N_G.v(t, theta0, delta0, r, mu, rho)
}

loss_death <- function(t, theta0, delta0, r, mu, rho){
  
  ifelse(t>0, integrate(numerator_delta, lower = 0, upper = t, t=t, theta0 = theta0, delta0 = delta0, r = r, mu = mu, rho = rho)$value/
     N_G.v(t, theta0, delta0, r, mu, rho),
     delta0)
}

loss_death_vec <- Vectorize(loss_death)


death_cells <- function(t, theta0, delta0, r, mu, rho){
  
  loss_death_vec(t, theta0, delta0, r, mu, rho) * N_G.v(t, theta0, delta0, r, mu, rho)
}


numerator_mu <-  function(a, t,theta0, delta0, r, mu, rho){
  
  mu * G_vec(a, t, theta0, delta0, r, mu, rho)/
    N_G.v(t, theta0, delta0, r, mu, rho)
}

loss_mu <- function(t, theta0, delta0, r, mu, rho){
  
  #integrate(numerator_mu, lower = 0, upper = t, t=t, theta0 = theta0, delta0 = delta0, r = r, mu = mu, rho = rho)$value/
    mu * N_G.v(t, theta0, delta0, r, mu, rho)
}

loss_mu_vec <-  Vectorize(loss_mu)

age_pred1 <- seq(0, 50)
age_pred2 <- seq(0, 150)
age_pred3 <- seq(0, 300)
age_pred4 <- seq(0, 450)

ts_pred <-  seq(0, 450)

ggplot()+
  geom_line(aes(ts_pred, N_G.v(ts_pred, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=1)+
  geom_line(aes(ts_pred, death_cells(ts_pred, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=2)+
  geom_line(aes(ts_pred, loss_mu_vec(ts_pred, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=4)+
  labs(title ="Counts of cells entering memory", y = NULL, x = "Host age(days)") + scale_y_log10()


ggplot()+
  geom_line(aes(age_pred1, normalised_death(age_pred1, 50, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=1) +
  geom_line(aes(age_pred2, normalised_death(age_pred2, 150, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=2) +
  geom_line(aes(age_pred3, normalised_death(age_pred3, 300, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=4) +
  geom_line(aes(age_pred4, normalised_death(age_pred4, 450, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=6) +
  labs(title ="Age dist of dying cells ", y = "Frequency", x = "Cell age(days)")
  
ggplot()+
  geom_line(aes(age_pred1, numerator_mu(age_pred1, 50, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=1) +
  geom_line(aes(age_pred2, numerator_mu(age_pred2, 150, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=2) +
  geom_line(aes(age_pred3, numerator_mu(age_pred3, 300, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=4) +
  geom_line(aes(age_pred4, numerator_mu(age_pred4, 450, 6, 0.05, 0.005, 0.01, 0.001)), size =1.5, col=6) +
  labs(title ="normalised cell age dist of cells entering memory", y = "Frequency", x = "Cell age(days)")


set.seed(1004)

ts_fit <- as.integer(c(runif(20, 5, 39),
            runif(50, 40, 400)))

counts_err <- c(rnorm(20, exp(13.5), exp(14.5)),
                rnorm(50, exp(13), exp(14.2)))

fake_data <- data.frame("host_age" = ts_fit, 
                        "counts" = N_G.v(ts_fit, 6, 0.08, 80) + counts_err) %>%
  arrange(host_age)


ts_pred <- seq(1, 400, 0.5)

ggplot()+
  geom_point(data = fake_data, aes(x = host_age, y = counts), size = 2) +
  geom_line(aes(ts_pred, N_G.v(ts_pred, 6, 0.08, 80)), size = 1.5) + scale_y_log10()
  
#parameters
v= v_cd4 #rate of thymic involution

#fitting log N_total & logit N_fd to the NCD4 data
ADT_NCD4_tx <- function(param, boot_data) { 
  theta0 <-param[1]
  l <- param[2]                     #parametrs to be estimated as part of a vector
  r <- param[3]            
  
  exp_WT.df <- boot_data%>%
    mutate(log_WT = log10(counts)) 
  
  pred_WT <- log10(N_G.v(exp_WT.df$host_age, theta0, l, r))
  
  ssqres1 <- sum((pred_WT - exp_WT.df$log_WT)^2)
  
  k  <- length(param)                #number of unknown parameters 
  n1 <- nrow(exp_WT.df)              #number of observations in dataset1
  
  #cost function
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multipltication of residual and transpose of residuals
  logl <-  - (n1/2)* log(ssqres1)  
  #stats for model validation
  AICcd4tx <<-  2*k - 2*logl
  
  return(-logl)     #since optim minimizes the function by default, ML
} 

ADT_N4Tx_fits <- optim(par=c(6, 0.03, 20), fn = ADT_NCD4_tx, boot_data = fake_data, hessian = TRUE)
ADT_N4Tx_fits
AICcd4tx


ggplot() +
  geom_point(data = fake_data, aes(host_age, y= counts), size = 2) +
  labs(title="Age-conditioning model", y="naive CD4 cell counts", x="mouse age (days)") + xlim(0, 400) +
  geom_line(aes(x= ts_pred, y=N_G.v(ts_pred, ADT_N4Tx_fits$par[1],  ADT_N4Tx_fits$par[2], ADT_N4Tx_fits$par[3])), color="#2548b4", size = 1.5)+
  scale_y_continuous(limits=c(1e6, 2e7),  trans="log10", breaks=c(1e7, 1e6, 1e8)) +
  theme_bw() + theme(axis.text = element_text(size = 32),
                     axis.title =  element_text(size = 22, face = "bold"),
                     panel.background = element_rect(colour = "black"),
                     plot.title = element_text(size=22, face = "bold",  hjust = 0.5),
                     legend.text = element_text(size=22, hjust = 0.5),
                     legend.title = element_blank())



params <- c("theta0" = 6, 'l' = 0.1, 'r' = 40)

ts_ode <- seq(1, 100)

foreach_time <- function(ts){
  
  ## ode function
  ode_func <-  function(s, state, parms){
    with(as.list(c(state, parms)),{
      #rate of change
      dX <-  -lambda_hill(s, l, r) * X
      
      #return the rate of change
      list(c(dX = theta_func(ts - s, theta0)))
      
    })  # end with(as.list ...
  }
  
  #initial conditions
  state <- c(X=0)
   
  # time points for which conc is reported
  # include the points where data is available
  times= c(0, ts)
  
  sol_ode <- ode(y= state, times=times, func = ode_func, parms = list(theta0 = 6, l = 0.01, r = 40))
  sol_ode[2, 2]
}

sol_df <- data.frame("time" = ts_ode)
for (i in ts_ode){
  
  sol_df[i,2] <- foreach_time(i)
}

int_df <- data.frame("time" = ts_ode,
                     "counts" = N_G.v(ts_ode, 6, 0.01, 40))

ggplot()+
  geom_line(data = sol_df, aes(time, V2), size =1.5, col =2)+
  geom_line(data = int_df, aes(time, counts), size =1.5, col =4)
  

foreach_time(100)























