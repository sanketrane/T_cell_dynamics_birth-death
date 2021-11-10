library(tidyverse)
library(deSolve)


# the function for total thymic output
theta_func <- function(Time, psi) {
  basl  = 2.5e6
  theta = 2e3
  n     = 2
  X     = 34
  q     = 3.5
  
  theta = basl + (theta * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q))))
  
  return(psi * theta)
}

### rate of turnover varying with cell age
delta_age <- function(age, delta0, r){
  
  delta0 * exp(-r * age)
}

lambda_age <- function(age, delta0, r, mu, rho){
  
  delta_age(age, delta0, r) + mu -  rho
}

# age distribution of cells
G_noabmt <- function(a, t, psi, delta0, r, mu, rho) {  
  
  theta_func(t-a, psi)  * exp(- integrate(lambda_age, lower = 0, upper = t, delta0 = delta0, r = r, mu = mu, rho = rho)$value)
}  

G_vec <-  Vectorize(G_noabmt)

N_G_noabmt <- function(t, psi, delta0, r, mu, rho){
  
  integrate(G_vec, lower = 0, upper = t, t = t, psi = psi, delta0 = delta0, r = r, mu = mu, rho = rho)$value
}

N_G.v <- Vectorize(N_G_noabmt)

numerator_delta <-  function(a, t,psi, delta0, r, mu, rho){
  
  (delta_age(a, delta0, r) * G_vec(a, t, psi, delta0, r, mu, rho))
}

normalised_death <- function(a, t,psi, delta0, r, mu, rho){
  
  (delta_age(a, delta0, r) * G_vec(a, t, psi, delta0, r, mu, rho))/
    N_G.v(t, psi, delta0, r, mu, rho)
}

loss_death <- function(t, psi, delta0, r, mu, rho){
  
  ifelse(t>0, integrate(numerator_delta, lower = 0, upper = t, t=t, psi = psi, delta0 = delta0, r = r, mu = mu, rho = rho)$value/
           N_G.v(t, psi, delta0, r, mu, rho),
         delta0)
}

loss_death_vec <- Vectorize(loss_death)


death_cells <- function(t, psi, delta0, r, mu, rho){
  
  loss_death_vec(t, psi, delta0, r, mu, rho) * N_G.v(t, psi, delta0, r, mu, rho)
}


numerator_mu <-  function(a, t,psi, delta0, r, mu, rho){
  
  mu * G_vec(a, t, psi, delta0, r, mu, rho)/
    N_G.v(t, psi, delta0, r, mu, rho)
}

loss_mu <- function(t, psi, delta0, r, mu, rho){
  
  #integrate(numerator_mu, lower = 0, upper = t, t=t, psi = psi, delta0 = delta0, r = r, mu = mu, rho = rho)$value/
  mu * N_G.v(t, psi, delta0, r, mu, rho)
}

loss_mu_vec <-  Vectorize(loss_mu)

age_pred1 <- seq(0, 50)
age_pred2 <- seq(0, 150)
age_pred3 <- seq(0, 300)
age_pred4 <- seq(0, 450)

ts_pred <-  seq(0, 450)

counts_df <- data.frame("host_age" = ts_pred,
                        "total" = N_G.v(ts_pred, 0.8, 0.05, 0.01, 0.005, 0.001),
                        "death" = death_cells(ts_pred, 0.8, 0.05, 0.01, 0.005, 0.001),
                        "differentation" = loss_mu_vec(ts_pred, 0.8, 0.05, 0.01, 0.005, 0.001)) %>%
  gather(-host_age, key = "process", value = "counts")

ggplot(counts_df)+
  geom_line(aes(host_age, counts, col = process), size =1.5)+
  labs(title ="Counts of cells", y = NULL, x = "Host age(days)") + scale_y_log10() 


ggplot()+
  geom_line(aes(age_pred1, normalised_death(age_pred1, 50, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=1) +
  geom_line(aes(age_pred2, normalised_death(age_pred2, 150, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=2) +
  geom_line(aes(age_pred3, normalised_death(age_pred3, 300, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=4) +
  geom_line(aes(age_pred4, normalised_death(age_pred4, 450, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=6) +
  labs(title ="Normalised Age dist of dying cells ", y = "Frequency", x = "Cell age(days)")  + scale_y_log10()

ggplot()+
  geom_line(aes(age_pred1, numerator_mu(age_pred1, 50, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=1) +
  geom_line(aes(age_pred2, numerator_mu(age_pred2, 150, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=2) +
  geom_line(aes(age_pred3, numerator_mu(age_pred3, 300, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=4) +
  geom_line(aes(age_pred4, numerator_mu(age_pred4, 450, 0.8, 0.03, 0.01, 0.01, 0.001)), size =1.5, col=6) +
  labs(title ="Normalised Age dist of cells entering memory", y = "Frequency", x = "Cell age(days)")+ scale_y_log10()

