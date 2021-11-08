library(rstan)
library(loo)
library(tidyverse)

Output_dir <- "output_csv"
m0 <- read_stan_csv(file.path(Output_dir, 'asmS49_deltavar_cd4.csv')) 
m1 <- read_stan_csv(file.path(Output_dir, 'asmS49_rhovar_cd4.csv')) 
m3 <- read_stan_csv(file.path(Output_dir, 'asmS52_deltavar_cd8.csv')) 
m7 <- read_stan_csv(file.path(Output_dir, 'asmS52_rhovar_cd8.csv')) 


## calculating an output from individual runs for the validation of a successful run
# loo-ic values
loo(extract_log_lik(m0, parameter_name = "log_lik_counts"), cores = 4)
loo(extract_log_lik(m1, parameter_name = "log_lik_counts"), cores = 4)
loo(extract_log_lik(m3, parameter_name = "log_lik_counts"), cores = 4)
loo(extract_log_lik(m7, parameter_name = "log_lik_counts"), cores = 4)

loo(extract_log_lik(m0, parameter_name = "log_lik_ki"), cores = 4)
loo(extract_log_lik(m1, parameter_name = "log_lik_ki"), cores = 4)
loo(extract_log_lik(m3, parameter_name = "log_lik_ki"), cores = 4)
loo(extract_log_lik(m7, parameter_name = "log_lik_ki"), cores = 4)


num_pars <- which(m0@model_pars == 'sigma_ki')
parameters_to_plot <- m0@model_pars[1: num_pars]

ptable <- monitor(as.array(m0, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
ptable1 <- monitor(as.array(m1, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table1 <- ptable1[1:num_pars, c(1, 3, 4, 8)]
ptable2 <- monitor(as.array(m2, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table2 <- ptable2[1:num_pars, c(1, 3, 4, 8)]
ptable3 <- monitor(as.array(m3, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table3 <- ptable3[1:num_pars, c(1, 3, 4, 8)]


##### main plots
Y0pred <- as.data.frame(m0, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

Y1pred <- as.data.frame(m1, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

Y3pred <- as.data.frame(m3, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

Y7pred <- as.data.frame(m7, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)


ki0pred <- as.data.frame(m0, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

ki1pred <- as.data.frame(m1, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

ki3pred <- as.data.frame(m3, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

ki7pred <- as.data.frame(m7, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

#### facet plot for individual age bins
ggplot() +
  geom_line(data = Y0pred, aes(x = timeseries, y = median), size =1.2, col = 4) +
  geom_line(data = Y1pred, aes(x = timeseries, y = median), size =1.2, col = 3) +
   #geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#43b1df", alpha = 0.15) +
  #geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#43b1df", alpha = 0.3) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_counts, col = dataset), alpha = 0.8, size = 2) +
  scale_color_manual(name=NULL, values=c(2, 6)) +
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8)) 

ggplot() +
  geom_line(data = ki0pred, aes(x = timeseries, y = median * 100), size =1.2, col = 4) +
  geom_line(data = ki1pred, aes(x = timeseries, y = median * 100), size =1.2, col = 3) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_kiprop * 100, col = dataset), alpha = 0.8, size = 2) +
  scale_color_manual(name=NULL, values=c(2, 6)) +
  labs(title= paste0('% of Ki67+ cell in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(0.3, 100), trans = 'log10') 



ggplot() +
  geom_line(data = Y3pred, aes(x = timeseries, y = median), size =1.2, col = 4) +
  geom_line(data = Y7pred, aes(x = timeseries, y = median), size =1.2, col = 3) +
  #geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#43b1df", alpha = 0.15) +
  #geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#43b1df", alpha = 0.3) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_counts, col = dataset), alpha = 0.8, size = 2) +
  scale_color_manual(name=NULL, values=c(2, 6)) +
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8)) 

ggplot() +
  geom_line(data = ki3pred, aes(x = timeseries, y = median * 100), size =1.2, col = 4) +
  geom_line(data = ki7pred, aes(x = timeseries, y = median * 100), size =1.2, col = 3) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_kiprop * 100, col = dataset), alpha = 0.8, size = 2) +
  scale_color_manual(name=NULL, values=c(2, 6)) +
  labs(title= paste0('% of Ki67+ cell in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(0.3, 100), trans = 'log10') 

