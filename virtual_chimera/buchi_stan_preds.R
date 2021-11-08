rm(list = ls()); gc()
library(rstan)
library(tidyverse)

theme_set(theme_bw())
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), 
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))
substrRight <- function(x, n){
  substr(x, 1, nchar(x)-n)
}

ageBMT_names <- c(`ageBMT_group1` = "7-9wks",
                  `ageBMT_group2` = "9-11wks", 
                  `ageBMT_group3` = "11-25wks")


## importing ASM stan model
ModelName <- 'neutralS49'
Population <- "cd4"
OutputDir <- file.path("output_csv")
SaveDir <- file.path('out_fit', Population, ModelName) 

stan_file <- paste0(substrRight(ModelName, 3), '_', Population, '.stan')
expose_stan_functions(file.path('stan_models/virtual_chimera', stan_file))


## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) 

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34))

data_counts <- chimera_data %>% 
  select(contains('age'), contains('counts'), contains('kiprop')) %>%
  mutate(dataset = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                          ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                 'ageBMT_group3'))) %>%
  bind_rows(ontogeny_data)


ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_counts, col = dataset), size =2.5) + 
  labs(x = 'Host age', y = NULL, title = paste0("Total counts of naive ", Population, " T cells")) +
  scale_y_log10(limits = c(1e5, 1e8)) 

ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_kiprop * 100, col = dataset), size =2.5) + 
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67+ cells in naive ", Population, " T cells")) +
  scale_y_log10(limits = c(0.3, 100)) 

ggplot() +
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = ki_dh_ratio, col = ageBMT_bin), size =2.5) + 
  labs(x = 'Host age', y = NULL, title = paste0("Ratio of donor to host Ki67+ cell in naive ", Population, " T cells")) +
  scale_y_log10(limits = c(0.3, 100)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + guides(col = F)



# vector of parameters
param_file <- file.path(SaveDir, paste0('params_', Population, '_', ModelName, '.csv'))
par_df <- data.frame(read.csv(param_file))[1:8, 'mean']
y1_0 <- par_df[5] * par_df[1] ### kappa0 * N0
y2_0 <- (1 - par_df[5]) * par_df[1] ### (1 - kappa0) * N0
parstan_vec <- c(par_df[4], par_df[2], par_df[3], par_df[6]) 

## init conditions at t0 = 5
init_cond_ont <- c(y1_0, y2_0)

## time course for predictions
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)


## stan solution onotgeny
ode_df <- solve_ode_ont(ts_pred_ont, init_cond_ont, parstan_vec)

stan_pred_df <- data.frame("age.at.S1K" = ts_pred_ont,
                            "y_pred" = matrix(unlist(ode_df), nrow = length(ode_df), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2,
         ki_prop = y_pred.1/(y_pred.1 + y_pred.2))



## stan solution chimeras
ta_counts <- solve_ode_ont(c(5, 40), init_cond_ont, parstan_vec)
init_cond_chi <- c(0, 0, ta_counts[[2]])

chivec1 <- sapply(ts_pred_chi1 - 54, Chi_spline)
chivec2 <- sapply(ts_pred_chi2 - 71, Chi_spline)
chivec3 <- sapply(ts_pred_chi3 - 97, Chi_spline)

ode_df1 <- solve_ode_chi(ts_pred_chi1, 54, init_cond_chi, parstan_vec)
ode_df2 <- solve_ode_chi(ts_pred_chi2, 71, init_cond_chi, parstan_vec)
ode_df3 <- solve_ode_chi(ts_pred_chi3, 97, init_cond_chi, parstan_vec)

stan_pred_df1 <- data.frame("age.at.S1K" = ts_pred_chi1,
                            "y_pred" = matrix(unlist(ode_df1), nrow = length(ode_df1), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         Nfd = (y_pred.1 + y_pred.2)/(total_counts * chivec1),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2),
         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
         ageBMT_bin = rep('ageBMT_group1', length(ts_pred_chi1)))

stan_pred_df2 <- data.frame("age.at.S1K" = ts_pred_chi2,
                            "y_pred" = matrix(unlist(ode_df2), nrow = length(ode_df2), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         Nfd = (y_pred.1 + y_pred.2)/ (total_counts * chivec2),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2),
         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
         ageBMT_bin = rep('ageBMT_group2', length(ts_pred_chi2)))

stan_pred_df3 <- data.frame("age.at.S1K" = ts_pred_chi3,
                            "y_pred" = matrix(unlist(ode_df3), nrow = length(ode_df3), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4,
         Nfd = (y_pred.1 + y_pred.2)/ (total_counts * chivec3),
         donor_ki = y_pred.1/(y_pred.1 + y_pred.2),
         host_ki = y_pred.3/(y_pred.3 + y_pred.4),
         ageBMT_bin = rep('ageBMT_group3', length(ts_pred_chi3)))


counts_pooled <- rbind(stan_pred_df1[c('age.at.S1K', 'ageBMT_bin', 'total_counts')],
                       stan_pred_df2[c('age.at.S1K', 'ageBMT_bin', 'total_counts')],
                       stan_pred_df3[c('age.at.S1K', 'ageBMT_bin', 'total_counts')])

Nfd_pooled <- rbind(stan_pred_df1[c('age.at.S1K', 'ageBMT_bin', 'Nfd')],
                       stan_pred_df2[c('age.at.S1K', 'ageBMT_bin', 'Nfd')],
                       stan_pred_df3[c('age.at.S1K', 'ageBMT_bin', 'Nfd')])

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file) 

kidonor_df <- data.frame("ageBMT_group1" = stan_pred_df1$donor_ki,
                         "ageBMT_group2" = stan_pred_df2$donor_ki,
                         "ageBMT_group3" = stan_pred_df3$donor_ki) %>%
  gather(key = ageBMT_bin, value = donor_ki) 

kihost_df <- data.frame("ageBMT_group1" = stan_pred_df1$host_ki, 
                        "ageBMT_group2" = stan_pred_df2$host_ki, 
                        "ageBMT_group3" = stan_pred_df3$host_ki) %>%
  gather(key = ageBMT_bin, value = host_ki)

kipred_df <- data.frame(kidonor_df, kihost_df) %>%
  select(-"ageBMT_bin.1") %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3)) %>%
  gather(-c(age.at.S1K, ageBMT_bin), key = subpop, value = ki_prop) 



pdf(file = file.path(SaveDir, "PredFig%03d.pdf"),
    width = 12, height = 9, onefile = FALSE, useDingbats = FALSE )

counts_plot <- ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_counts, col = dataset), size =2.5) + 
  geom_line(data = stan_pred_df, aes(x = age.at.S1K, y = total_counts), size =1.25, col = 6) + 
  geom_line(data = counts_pooled, aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size =1.25) + 
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3, 6), labels = ageBMT_names)+
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  #scale_x_log10(limits = c(5, 450), breaks = c(10, 30, 100, 300, 600)) +
  myTheme  + guides(col =F)

Nfd_plot <- ggplot() +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size=1.3) +
  geom_line(data = Nfd_pooled, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size =1.25) + 
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size =2.5) + 
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3), labels = ageBMT_names)+
  labs(x = 'Host age', y = NULL, title = paste0("Normalised donor fraction in naive ", toupper(Population), " T cells")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  #scale_x_log10(limits = c(50, 450), breaks = c(10, 30, 100, 300, 450)) +
  myTheme + theme(legend.position = c(0.85, 0.15))


ki_facets <- ggplot() +
  geom_line(data = kipred_df, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size =1.25) +
  geom_point(data = ki_data, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size =2) + 
  scale_color_manual(name = NULL, values = c(2, 4), labels = c('Donor', "Host")) +
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67high cells in naive ", toupper(Population), " subset")) +
  scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100)) +
  #scale_x_log10(limits = c(5, 450), breaks = c(10, 30, 100, 300, 600)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + theme(legend.position = c(0.95, 0.9)) 

toprow <- cowplot::plot_grid(counts_plot,Nfd_plot,  labels = c("A", "B"), ncol = 2)
cowplot::plot_grid(toprow, ki_facets,  labels = c("", "C"), nrow = 2)



dev.off()


### R sol
library(deSolve)
theta_R <- function(t, parms){
  psi = parms[1]
  
  psi * sp_numbers(t)
}

chi_R <- function(t){
  chiEst = 0.847543332;
  qEst = 0.050944623;
  if (t < 0){chi = 0} else {
    chi = chiEst * (1 - exp(-qEst * t));
  }
  return(chi)
}

eps_R <- function(t){
  eps_0 = 0.14965320; eps_f = 0.03470231; A = 3.43078629;
  exp(-eps_f * (t + A)) + eps_0
}

plot(eps_R(ts_pred_chi1)~ts_pred_chi1)

beta = 1/3.5

## ode function
ode_func <-  function(t, state, parms){
  with(as.list(c(state, parms)),{
    
    # ki hi donor
    dy1 <- theta_R(t, parms) * chi_R(t - tBMT) * eps_R(t) + rho * (2 * y2 + y1) -
      (beta + delta) * y1
    
    # ki lo donor
    dy2 <- theta_R(t, parms) * chi_R(t - tBMT) * (1 - eps_R(t)) + beta * y1 - 
      (rho + delta) * y2
    
    # ki hi host
    dy3 <- theta_R(t, parms) * (1 - chi_R(t - tBMT)) * eps_R(t) + rho * (2 * y4 + y3) -
      (beta + delta) * y3
    
    # ki lo host
    dy4 <- theta_R(t, parms) * (1 - chi_R(t - tBMT)) * (1 - eps_R(t)) + beta * y3 - 
      (rho + delta) * y4
    
    #return the rate of change
    list(c(dy1, dy2, dy3, dy4))
    
  })  # end with(as.list ...
}

#initial conditions
state <- init_cond
names(state) <- c("y1", "y2", "y3", "y4")

parms_vec <- c('psi' = 0.356223345, 'delta' = 0.019323868, 'rho' = 0.001937711)
tBMT = 58
sol_ode1 <- ode(y= state, times=ts_pred_chi1, func = ode_func, parms = parms_vec)
tBMT = 75
sol_ode2 <- ode(y= state, times=ts_pred_chi2, func = ode_func, parms = parms_vec)
tBMT = 105
sol_ode3 <- ode(y= state, times=ts_pred_chi3, func = ode_func, parms = parms_vec)

chivec1 <- chi_R(ts_pred_chi1 - 58)
chivec2 <- chi_R(ts_pred_chi2 - 75)
chivec3 <- chi_R(ts_pred_chi3 - 105)

sol_df1 <-  data.frame(sol_ode1) %>%
  mutate(age.at.S1K = ts_pred_chi1,
         total_counts = y1 + y2 + y3 + y4,
         donor_ki = y1/ (y1 + y2),
         host_ki = y3/(y3 + y4),
         fd = (y1 + y2)/(total_counts * chivec1),
         AgeBMTbin = rep('age_group1', length(ts_pred_chi1)))#%>% na.omit

sol_df2 <-  data.frame(sol_ode2) %>%
  mutate(age.at.S1K = ts_pred_chi2,
         total_counts = y1 + y2 + y3 + y4,
         donor_ki = y1/ (y1 + y2),
         host_ki = y3/(y3 + y4),
         fd = (y1 + y2)/(total_counts * chivec2),
         AgeBMTbin = rep('age_group2', length(ts_pred_chi2))) #%>% na.omit

sol_df3 <-  data.frame(sol_ode3) %>%
  mutate(age.at.S1K = ts_pred_chi3,
         total_counts = y1 + y2 + y3 + y4,
         donor_ki = y1/ (y1 + y2),
         host_ki = y3/(y3 + y4),
         fd = (y1 + y2)/(total_counts * chivec3),
         AgeBMTbin = rep('age_group3', length(ts_pred_chi3))) #%>% na.omit


counts_R <- rbind(sol_df1[c('age.at.S1K', 'AgeBMTbin', 'total_counts')],
                  sol_df2[c('age.at.S1K', 'AgeBMTbin', 'total_counts')],
                  sol_df3[c('age.at.S1K', 'AgeBMTbin', 'total_counts')])


ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_counts, col = dataset), size =2.5) + 
  #geom_line(aes(x = ts_pred_ont, y = ntotal_vec), size =1.25, col = 5) + 
  geom_line(data = counts_R, aes(x = age.at.S1K, y = total_counts, col = AgeBMTbin),
            size =1.25) + 
  labs(x = 'Host age', y = NULL, title = paste0("Total counts of naive ", Population, " T cells")) +
  scale_y_log10(limits = c(1e5, 1e8)) + 
  #scale_x_log10(limits = c(5, 450), breaks = c(10, 30, 100, 300, 600)) +
  myTheme


