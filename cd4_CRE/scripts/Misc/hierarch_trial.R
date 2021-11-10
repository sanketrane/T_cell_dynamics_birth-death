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
### data for fitting
naive_df <- tstmp_df %>%
  select('id', 'age_anim', 'age_group', 'naive_labelled')  %>%
  rename(mouse_id = id, day_stamp = age_group) %>%
  mutate(group_id = ifelse(day_stamp == 1, 1,
                           ifelse(day_stamp == 7, 2,
                                  ifelse(day_stamp == 28, 3,
                                         ifelse(day_stamp == 56, 4, 5)))),
         t0_group = ifelse(day_stamp == 1, 14,
                           ifelse(day_stamp == 7, 28,
                                  ifelse(day_stamp == 28, 42,
                                         ifelse(day_stamp == 56, 70, 189))))) %>%
  #filter(group_id == c(1, 2)) %>%
  arrange(mouse_id)


ggplot()+
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                 col = as.factor(mouse_id)), size =3) + guides(col=F)+
  #geom_line(data= naive_df,  aes(x = age_anim, y = naive_labelled,
   #                               col = as.factor(mouse_id)), size =1.2) +
  scale_color_discrete(name = "Mouse ID", guide = guide_legend(nrow = 4)) +
  scale_y_log10(limits = c(5e3, 1e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 120), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_rect(fill = NA, colour = NA))


write.csv(naive_df, file = file.path("data", paste0("naive_labelled.csv")))



#library(rstan)
library(cmdstanr)
library(parallel)
# stan inputs -------------------------------
### data 
data_list <- list(
  K = length(unique(naive_df$mouse_id)),
  J = length(unique(naive_df$group_id)),
  N = nrow(naive_df),
  age_anim = naive_df$age_anim,
  naive_labelled = naive_df$naive_labelled,
  mouse_id = naive_df$mouse_id,
  group_id = naive_df$group_id,
  day_stamp = unique(naive_df$day_stamp),
  t0_group = unique(naive_df$t0_group),
  numPred = 101,
  ts_pred = seq(0, 100)
)

init_list <- function(){
  list(
    mu_lambda = exp(rnorm(1, log(0.05), 0.1)),
    mu_N0     = rnorm(1, 5.5, 0.1),
    r_lamda   = exp(rnorm(1, log(0.01), 0.1)),
    sigma_counts = runif(1, 0.1, 0.3),
    sigma_N0   = runif(1, 0, 0.5),
    sigma_lambda    = runif(1, 0, 0.01)
  )
}

parameters_to_plot <- c("r_lambda",
                        "N0",
                        "mu_N0",
                        "sigma_N0",
                        "lambda",
                        "mu_lambda",
                        "sigma_lambda",
                        "sigma_counts"
)

other_rvs <- c("log_lik_counts")

### params  and interesting quantities to track in the model
parameters <- c(parameters_to_plot, other_rvs)

### compile model
ModelName <- "hierarch_pure_ASM"
model_file <- file.path("stan_models", paste0(ModelName, '.stan'))
mod <- cmdstan_model(model_file)

### sampling params to fit model to data
fit <- mod$sample(data = data_list, init = init_list, 
                  iter_warmup = 500, iter_sampling = 1000,
                  chains = 4, parallel_chains = 4,
                  save_warmup = T, refresh = 100, 
                  adapt_delta = 0.9)

# saving the stan fit object for individual runs as 
stanfit_save <- rstan::read_stan_csv(fit$output_files())
output_filename <- file.path("Fits_save",  paste0(ModelName, ".rds"))
write_rds(stanfit_save, file = file.path(output_filename))

#pairs(stanfit_save, pars = parameters_to_plot)


### parameters table

ptable <- monitor(as.array(stanfit_save, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:37, c(1, 3, 4, 8)]
out_table

# loo-ic values
loo_loglik <- loo::extract_log_lik(stanfit_save, parameter_name = "log_lik_counts", merge_chains = TRUE)
loo_ic <- loo(loo_loglik,  save_psis = FALSE, cores = 4)
loo_ic$estimates

stanfit_save <- read_rds(file.path("Fits_save",  paste0(ModelName, ".rds")))

lambda_pred <- as.data.frame(stanfit_save, pars = "lambda") %>%
  'colnames<-' (seq(1, 66)) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) 


ggplot()+
  geom_boxplot(data = lambda_pred, aes(x=key, y=value, fill=key)) +
  scale_fill_manual(values = rep('#6a9ea4', 66)) + 
  guides(fill = FALSE) +
  labs(title = 'Variation in loss rate within mice', x = 'mouse ID', y= 'Loss rate of cells of age 0 ')


N0_pred <- as.data.frame(stanfit_save, pars = "N0") %>%
  'colnames<-' (seq(1, 5)) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) 


ggplot()+
  geom_boxplot(data = N0_pred, aes(x=key, y=10^value, fill=key)) +
  scale_fill_manual(values = rep('#6a9ea4', 5)) + 
  guides(fill = FALSE) + scale_y_log10(limits = c(1e5, 2e6))+
  labs(title = 'Variation in N0 within cohorts', x = 'Group ID', y= 'N0')


# vector of parameters
parstan1_vec <- c(mean(out_table$mean[2:13]), out_table$mean[70], out_table$mean[1])
parstan2_vec <- c(mean(out_table$mean[14:30]), out_table$mean[70], out_table$mean[1])
parstan3_vec <- c(mean(out_table$mean[31:42]), out_table$mean[70], out_table$mean[1])
parstan4_vec <- c(mean(out_table$mean[43:54]), out_table$mean[70], out_table$mean[1])
parstan5_vec <- c(mean(out_table$mean[55:67]), out_table$mean[70], out_table$mean[1])
rdata1_vec <- c(1, 14)
rdata2_vec <- c(7, 28)
rdata3_vec <- c(28, 42)
rdata4_vec <- c(56, 70)
rdata5_vec <- c(175, 189)

## stan solution
sim1 <- sapply(tseq1,  N_total_time, parms = c(parstan1_vec, rdata1_vec))
sim2 <- sapply(tseq2,  N_total_time, parms = c(parstan2_vec, rdata2_vec))
sim3 <- sapply(tseq3,  N_total_time, parms = c(parstan3_vec, rdata3_vec))
sim4 <- sapply(tseq4,  N_total_time, parms = c(parstan4_vec, rdata4_vec))
sim5 <- sapply(tseq5,  N_total_time, parms = c(parstan5_vec, rdata5_vec))

## labels for individual age cohorts
age_cohorts <- c(rep('Cohort-1', 100), rep('Cohort-2', 100), rep('Cohort-3', 100), 
                 rep('Cohort-4', 100), rep('Cohort-5', 100))
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


#ggsave(paste0(ModelName, '_', Population, 'simplot2.pdf'), width = 6, height = 4.5, units = "in")













