## clearing the environment
rm(list = ls())  
gc()    

## loading libraries
library(rstan)
library(loo)
library(tidyverse)

myTheme <- theme(text = element_text(size = 18), axis.text = element_text(size = 12),
                 axis.title =  element_text(size = 18, face = "bold"),
                 plot.title = element_text(size=18,  hjust = 0.5, face = "bold"))

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

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

log10minorbreaks = as.numeric(1:10 %o% 10^(4:8))


## Read model fit files
PureASM_N0grp_L0mouse_fit <- read_rds(file.path("Fits_save",  paste0("PureASM_N0grp_L0mouse", ".rds")))
PureASM_N0mouse_L0grp_fit <- read_rds(file.path("Fits_save",  paste0("PureASM_N0mouse_L0grp", ".rds")))
PureASM_N0mouse_L0mouse_fit <- read_rds(file.path("Fits_save",  paste0("PureASM_N0mouse_L0mouse", ".rds")))
PureASM_N0mouse_fit <- read_rds(file.path("Fits_save",  paste0("PureASM_N0mouse", ".rds")))
APHA_N0mouse_Sigmoid2_fit <- read_rds(file.path("Fits_save",  paste0("APHA_N0mouse_Sigmoid2", ".rds")))
APHA_N0mouse_Sigmoid5_fit <- read_rds(file.path("Fits_save",  paste0("APHA_N0mouse_Sigmoid5", ".rds")))
DDM_N0mouse_fit <- read_rds(file.path("Fits_save", 
                                      paste0("DDM_N0mouse", ".rds")))
DDM_N0mouse_L0grp_fit <- read_rds(file.path("Fits_save", 
                                      paste0("DDM_N0mouse_L0cohort", ".rds")))
DDM_N0mouse_L0mouse_fit <- read_rds(file.path("Fits_save", 
                                      paste0("DDM_N0mouse_L0mouse", ".rds")))

### LOOIC estimates
Loo_PureASM_N0grp_L0mouse <- loo(extract_log_lik(PureASM_N0grp_L0mouse_fit,
                                parameter_name = "log_lik_counts", merge_chains = TRUE),
                                 save_psis = FALSE, cores = 4)
Loo_PureASM_N0mouse_L0grp <- loo(extract_log_lik(PureASM_N0mouse_L0grp_fit,
                                             parameter_name = "log_lik_counts", merge_chains = TRUE),
                                 save_psis = FALSE, cores = 4)
Loo_PureASM_N0mouse_L0mouse <- loo(extract_log_lik(PureASM_N0mouse_L0mouse_fit,
                                             parameter_name = "log_lik_counts", merge_chains = TRUE),
                                   save_psis = FALSE, cores = 4)
Loo_PureASM_N0mouse <- loo(extract_log_lik(PureASM_N0mouse_fit,
                                             parameter_name = "log_lik_counts", merge_chains = TRUE),
                           save_psis = FALSE, cores = 4)

Loo_APHA_N0mouse_Sigmoid2 <- loo(extract_log_lik(APHA_N0mouse_Sigmoid2_fit,
                                                 parameter_name = "log_lik_counts", merge_chains = TRUE),
                                 save_psis = FALSE, cores = 4)

Loo_APHA_N0mouse_Sigmoid5 <- loo(extract_log_lik(APHA_N0mouse_Sigmoid5_fit,
                                                 parameter_name = "log_lik_counts", merge_chains = TRUE),
                                 save_psis = FALSE, cores = 4)


Loo_DDM_N0mouse <- loo(extract_log_lik(DDM_N0mouse_fit,
                                                 parameter_name = "log_lik_counts", merge_chains = TRUE),
                                 save_psis = FALSE, cores = 4)


Loo_DDM_N0mouse_L0grp <- loo(extract_log_lik(DDM_N0mouse_L0grp_fit,
                                       parameter_name = "log_lik_counts", merge_chains = TRUE),
                       save_psis = FALSE, cores = 4)

Loo_DDM_N0mouse_L0mouse <- loo(extract_log_lik(DDM_N0mouse_L0mouse_fit,
                                       parameter_name = "log_lik_counts", merge_chains = TRUE),
                       save_psis = FALSE, cores = 4)


model_comparison <- loo_compare(Loo_PureASM_N0grp_L0mouse, Loo_PureASM_N0mouse_L0grp,
             Loo_PureASM_N0mouse_L0mouse, Loo_PureASM_N0mouse,
             Loo_APHA_N0mouse_Sigmoid2, Loo_APHA_N0mouse_Sigmoid5)

print(model_comparison, simplify = FALSE, digits = 2)

model_list <- list("ASM_N0mouse" = Loo_PureASM_N0mouse,
                   "ASM_N0mouse_L0mouse" = Loo_PureASM_N0mouse_L0mouse,
                   "ASM_N0mouse_L0grp" = Loo_PureASM_N0mouse_L0grp,
                   "ASM_N0mouse_L0hostage" = Loo_APHA_N0mouse_Sigmoid5,
                   "DDM_N0mouse" = Loo_DDM_N0mouse,
                   "DDM_N0mouse_L0mouse" = Loo_DDM_N0mouse_L0mouse,
                   "DDM_N0mouse_L0grp" = Loo_DDM_N0mouse_L0grp
                   )
model_comparison_new <- loo_compare(model_list)

row.names(model_comp)

model_wts <- loo_model_weights(model_list, method = 'pseudobma')
model_comp <- data.frame(loo_compare(model_list)[,c(1, 2 , 7, 5)]) %>%
  mutate(delta_looic = looic - min(looic))

require(formattable)
export_table <- formattable(model_comp)

model_comp

model_compare <- function(tabl){
  # delta loo-ic
  deltaloo_list <- tabl$looic - min(tabl$looic)
  ## akaike wt
  akaikewt_numerator <- function(x) exp(-0.5 * x)
  akaikewt_list <- sapply(deltaloo_list, akaikewt_numerator) * 100/sum(exp(- 0.5 * deltaloo_list))
  
  export_table <- data.frame('deltaelpd' = round(tabl$elpd_diff, 2),
                             'deltaSE' = round(tabl$se_diff, 2),
                             'deltaloo' = round(deltaloo_list, 2),
                             'p_loo' = round(tabl$p_loo, 2),
                             'Akaike_wt' = round(akaikewt_list, 2))
  
  colnames(export_table)[1:5] <-  c(paste0('\u0394', 'elpd'), 
                                    paste0('\u0394', 'SE_elpd'), 
                                    paste0('\u0394', 'LooIC'), 
                                    paste0('Effective params'), 
                                    paste0('Akaike weight ', '\u0025'))
  
  return(export_table)
}

looic_table <- model_compare(model_comp)
export_table <- formattable(looic_table)
model_compare <- function(looiclist1){
  # delta loo-ic
  deltaloo_list <- looiclist1 - min(looiclist1)
  ## akaike wt
  akaikewt_numerator <- function(x) exp(-0.5 * x)
  akaikewt_list <- sapply(deltaloo_list, akaikewt_numerator) * 100/sum(exp(- 0.5 * deltaloo_list))
  
  export_table <- data.frame('deltaloo' = round(deltaloo_list, 2),
                             'Akaike_wt' = round(akaikewt_list, 2))
  colnames(export_table)[1:2] <-  c(paste0('\u0394', 'LooIC'),  paste0('Akaike weight ', '\u0025'))
  
  return(export_table)
}


looic_table <- model_compare(model_comp)

require(formattable)
export_table <- formattable(looic_table,
                            align = c('c', 'c', 'c'),
                            list(`p_loo` =
                                   formatter("span", style = ~ style(font.size='16px')),
                                 `elpd_diff` =
                                   formatter("span", style = x ~ ifelse(x == 0.00,
                                                                        style(color = customGreen, font.size='16px'),
                                                                        ifelse(x >=- 5.0, style(color = customOrange, font.size='16px'), 
                                                                               ifelse(x <100, style(font.size='16px'), NA)))))
)

export_table

#### predictions
model_fit <- PureASM_N0mouse_L0grp_fit#APHA_N0mouse_Sigmoid5_fit#

### parameters table
num_pars <- 6 + 66 

parameters_to_plot <- model_fit@model_pars[1:8]

ptable <- monitor(as.array(model_fit, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
View(out_table)


#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv") 
#model_fit <- stanfit_save
## prediction points
tseq1 <- seq(14, 120, length.out = 100)
tseq2 <- seq(28, 120, length.out = 100)
tseq3 <- seq(42, 100, length.out = 100)
tseq4 <- seq(70, 160, length.out = 100)
tseq5 <- seq(189, 280, length.out = 100)

##### main plots
Y1pred <- as.data.frame(model_fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = tseq1,
            'group_id' = rep(1, 100))

##### main plots
Y2pred <- as.data.frame(model_fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = tseq2,
            'group_id' = rep(2, 100))

##### main plots
Y3pred <- as.data.frame(model_fit, pars = "y3_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = tseq3,
            'group_id' = rep(3, 100))

##### main plots
Y4pred <- as.data.frame(model_fit, pars = "y4_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = tseq4,
            'group_id' = rep(4, 100))


##### main plots
Y5pred <- as.data.frame(model_fit, pars = "y5_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = tseq5,
            'group_id' = rep(5, 100))

naive_sim_df <- rbind(Y1pred, Y2pred, Y3pred, Y4pred, Y5pred)
naive_df$group_id <- as.factor(naive_df$group_id)
levels(naive_df$group_id) <- seq(1, 5)#c('Cohort-1', "Cohort-2", 'Cohort-3', "Cohort-4",'Cohort-5')


ggplot()+
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(group_id)), alpha=0.8, size =3) + 
  geom_line(data = naive_sim_df, aes(x = timeseries, y = median,
                                     col = as.factor(group_id)), size =1.5) + 
  geom_ribbon(data = naive_sim_df, aes(x = timeseries, ymin = lb, ymax = ub,
                                       fill = as.factor(group_id)), alpha = 0.4) +
  scale_color_discrete(name = 'Cohort', guide = guide_legend(nrow = 5)) +guides(fill = F) +
  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
  myTheme +
 # facet_wrap(~ as.factor(group_id), scales = 'free_x') + guides(col = F) +
  theme(#legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = NA, colour = NA))


#ggplot()+
#  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
#                                  col = as.factor(group_id)), size =2.2) + 
#  geom_line(data = naive_sim_df, aes(x = timeseries, y = median),
#                                     col = 'navy', size =1.5) + 
#  geom_ribbon(data = naive_sim_df, aes(x = timeseries, ymin = lb, ymax = ub),
#                                       fill = 'blue', alpha = 0.4) +
#  scale_color_brewer(palette = 2)+
#  scale_color_discrete(name = NULL, guide = guide_legend(nrow = 3)) +guides(fill = F) +
#  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
#  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
#  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
#  myTheme +
#  facet_wrap(~ as.factor(group_id), scales = 'free_x') + guides(col = F) +
#  theme(#legend.position = c(0.85, 0.87),
#        legend.background = element_rect(fill = NA, colour = NA))



ts_pred <- seq(1, 200, length.out = 200)
lambda0_hostage_pred <- as.data.frame(model_fit, pars = "lambda0_hostage_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

ggplot()+
  geom_line(data = lambda0_hostage_pred, aes(x = timeseries, y = median), size =1.5) + 
  geom_ribbon(data = lambda0_hostage_pred, aes(x = timeseries, ymin = lb, ymax = ub), alpha = 0.4)+
  scale_x_continuous(limits = c(0.7, 240), trans="log10", breaks = c(1, 3, 10, 30, 100, 300)) +
  guides(fill = FALSE) + myTheme +
  labs(title = bquote('Variation in ' ~ lambda(0) ~  ' with host age'),
       x = 'Host age', y= bquote( ~ lambda(0) ~ 'value'))


#lambda0_cd8ASMDeltavar <- c(0.129770922)


## mean host age when cohorts enetered the periphery (initial age dist => 6days)  
L0_meanAge <- unique(naive_df$day_stamp) #+ 3
lambda_pred <- as.data.frame(model_fit, pars = "lambda") %>%
  'colnames<-' (seq(1, 5)) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  mutate(group_id = ifelse(as.numeric(key) ==1, L0_meanAge[1],
                           ifelse(as.numeric(key) == 2, L0_meanAge[2],
                                  ifelse(as.numeric(key) == 3, L0_meanAge[3],
                                         ifelse(as.numeric(key) == 4, L0_meanAge[4], L0_meanAge[5])))))

lambda_pred22<- as.data.frame(model_fit, pars = "lambda") %>%
  'colnames<-' (seq(1, 5)) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  rename(cohort = key) %>%
  mutate(Host_age = ifelse(as.numeric(cohort) ==1, L0_meanAge[1],
                           ifelse(as.numeric(cohort) == 2, L0_meanAge[2],
                                  ifelse(as.numeric(cohort) == 3, L0_meanAge[3],
                                         ifelse(as.numeric(cohort) == 4, L0_meanAge[4], L0_meanAge[5])))))


#write.csv(lambda_pred22, 'lambda_var_cohorts.csv', row.names = F)

ggplot()+
  geom_violin(data = lambda_pred, aes(x=key, y=value, fill= as.factor(group_id))) +
  guides(fill = FALSE) + myTheme +
  labs(title = bquote('Variation in ' ~ lambda ~  ' within cohorts'),
       x = 'Cohort ID', y= bquote( ~ lambda(a == 0)))


ggplot()+
  geom_boxplot(data = lambda_pred, aes(x=group_id, y=value, fill= as.factor(key))) +
  geom_line(data = lambda0_hostage_pred, aes(x = timeseries, y = median), col='#a500ff', size =1.5) + 
  geom_ribbon(data = lambda0_hostage_pred, aes(x = timeseries, ymin = lb, ymax = ub),
              fill='#a500ff',alpha = 0.2)+
  scale_x_continuous(limits = c(0.7, 240), trans="log10", breaks = c(1, 3, 10, 30, 100, 300)) +
  scale_fill_discrete(name = 'Cohort') +
  #guides(fill = FALSE) + 
  myTheme +
  labs(title = bquote('Variation in ' ~ lambda ~  ' with host age'),
       x = 'Host age', y= bquote( lambda(t, a==0)))


ggplot()+
  geom_boxplot(data = lambda_pred, aes(x=group_id, y=value, fill= as.factor(group_id))) +
  #geom_line(aes(x=ts_pred, y=lambda_vec),  alpha=0.8, col = 'blue', size = 2)+
  #geom_line(aes(x=ts_pred, y=lambda_vec2), alpha=0.8, col = 'blue', size = 2)+
  #geom_line(aes(x=ts_pred, y=lambda_vec4), alpha=0.8, col = 6, size = 2)+
  scale_x_continuous(limits = c(0.7, 240), trans="log10", breaks = c(1, 3, 10, 30, 100, 300)) +
  guides(fill = FALSE) + myTheme +
  labs(title = bquote('Variation in ' ~ lambda ~  ' within cohorts'),
       x = 'Host age', y= bquote( ~ lambda(a == 0)))



lambda_form <- function(Time, l0, q){
  l0 * (1 + (1/(1+ (Time/q)^3)))
}


lambda_form2 <- function(Time, l0, q){
  
  l0 * exp(-q * Time)
}

lambda_form3 <- function(Time, l0, q1, q2){
  l0 * (1 + (q2/(1+ (Time/q1)^3)))
}

lambda_form4 <- function(Time, l0, q1, q2){
  l0 * (1 + (q2/(1+ (Time/q1)^4)))
}

ts_pred <- 10^ seq(0, 2.3, length.out = 100)
lambda_form5 <- function(Time, l0, q1, q2){
  (l0*(1 + q2/(1  + (log10(Time)/q1)^15))) #* (1+ eps - eps*exp(-5*log10(Time)))
  }

lambda_vec5 <- lambda_form5(ts_pred, 0.026, 1.51086, 1.48841)


ggplot()+
  geom_boxplot(data = lambda_pred, aes(x=group_id, y=value, fill= as.factor(group_id))) +
  #geom_line(aes(x=ts_pred, y=lambda_vec),  alpha=0.8, col = 'blue', size = 2)+
  #geom_line(aes(x=ts_pred, y=lambda_vec2), alpha=0.8, col = 'blue', size = 2)+
  geom_line(aes(x=ts_pred, y=lambda_vec5), alpha=0.8, col = 'blue', size = 2)+
  scale_x_continuous(limits = c(0.7, 240), trans="log10", breaks = c(1, 3, 10, 30, 100, 300)) +
  guides(fill = FALSE) + myTheme +
  labs(title = bquote('Variation in ' ~ lambda(0) ~  ' within cohorts'),
       x = 'Host age', y= bquote( ~ lambda(t, a=0) ~ 'value'))


lambda_vec <- lambda_form(ts_pred, 0.03082035, 31.14837760)
lambda_vec2 <- lambda_form2(ts_pred, 0.060489689, 0.006658414)

nls_fit <- nls(value ~ lambda_form5(group_id, l0, q1, q2, eps), data = lambda_pred,
               start = list('l0' = 0.06/2.5, 'q1' = 1.5, 'q2' = 1.5, 'eps' = 0.1))


nls_fit


N0_pred <- as.data.frame(model_fit, pars = "N0") %>%
  'colnames<-' (seq(1, 66)) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  mutate(group_id = ifelse(as.numeric(key) <=12, 'Cohort-1',
                           ifelse(as.numeric(key) <= 29, 'Cohort-2',
                                  ifelse(as.numeric(key) <=41, 'Cohort-3',
                                         ifelse(as.numeric(key) <= 53, 'Cohort-4', 'Cohort-5')))))


ggplot()+
  geom_boxplot(data = N0_pred, aes(x=key, y=10^value, fill=group_id)) +
  #scale_fill_manual(values = rep('#4890ab', 66)) + 
  guides(fill = FALSE) +  myTheme +
  scale_y_continuous(trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(title = bquote('Variation in ' ~ N(0) ~  ' within mice'),
       x = 'Mouse ID', y= bquote( ~ N(0) ~ 'count'))









