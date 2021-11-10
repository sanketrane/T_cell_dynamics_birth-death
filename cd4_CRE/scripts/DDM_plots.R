## clearing the environment
rm(list = ls())  
gc()    

## loading libraries
library(rstan)
library(loo)
library(tidyverse)

myTheme <- theme(text = element_text(size = 15), axis.text = element_text(size = 14),
                 axis.title =  element_text(size = 15, face = "bold"),
                 plot.title = element_text(size=15,  hjust = 0.5, face = "bold"))

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
DDM_N0mouse_fit <- read_rds(file.path("Fits_save", 
                                      paste0("DDM_N0mouse_L0cohort", ".rds")))

#### predictions
model_fit <- DDM_N0mouse_fit

### parameters table
num_pars <- 74 

parameters_to_plot <- model_fit@model_pars[1:8]

ptable <- monitor(as.array(model_fit, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
View(out_table)


#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv") 
model_fit <- stanfit_save
## prediction points
tseq1 <- seq(14, 120, length.out = 200)
tseq2 <- seq(28, 120, length.out = 200)
tseq3 <- seq(42, 100, length.out = 200)
tseq4 <- seq(70, 160, length.out = 200)
tseq5 <- seq(189, 280, length.out = 200)

##### main plots
Y1pred <- as.data.frame(model_fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("age_anim" = tseq1,
            'group_id' = rep(1, 200))

##### main plots
Y2pred <- as.data.frame(model_fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("age_anim" = tseq2,
            'group_id' = rep(2, 200))

##### main plots
Y3pred <- as.data.frame(model_fit, pars = "y3_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("age_anim" = tseq3,
            'group_id' = rep(3, 200))

##### main plots
Y4pred <- as.data.frame(model_fit, pars = "y4_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("age_anim" = tseq4,
            'group_id' = rep(4, 200))


##### main plots
Y5pred <- as.data.frame(model_fit, pars = "y5_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("age_anim" = tseq5,
            'group_id' = rep(5, 200))

naive_sim_df <- rbind(Y1pred, Y2pred, Y3pred, Y4pred, Y5pred)
naive_df$group_id <- as.factor(naive_df$group_id)
levels(naive_df$group_id) <- seq(1, 5)#c('Cohort-1', "Cohort-2", 'Cohort-3', "Cohort-4",'Cohort-5')


ggplot()+
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(group_id)), alpha=0.8, size =3) + 
  geom_line(data = naive_sim_df, aes(x = age_anim, y = median,
                                     col = as.factor(group_id)), size =1.5) + 
  geom_ribbon(data = naive_sim_df, aes(x = age_anim, ymin = lb, ymax = ub,
                                       fill = as.factor(group_id)), alpha = 0.4) +
  scale_color_discrete(name = 'Cohort', guide = guide_legend(nrow = 5)) +guides(fill = F) +
  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
  myTheme +
  # facet_wrap(~ as.factor(group_id), scales = 'free_x') + guides(col = F) +
  theme(#legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = NA, colour = NA))


## mean host age when cohorts enetered the periphery (initial age dist => 6days)  
L0_meanAge <- unique(naive_df$day_stamp) #+ 3
lambda_pred <- as.data.frame(model_fit, pars = "lambda") %>%
  'colnames<-' (seq(1, 66)) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  mutate(group_id = ifelse(as.numeric(key) <=12, 'Cohort-1',
                           ifelse(as.numeric(key) <= 29, 'Cohort-2',
                                  ifelse(as.numeric(key) <=41, 'Cohort-3',
                                         ifelse(as.numeric(key) <= 53, 'Cohort-4', 'Cohort-5')))))#%>%
  #mutate(group_id = ifelse(as.numeric(key) ==1, L0_meanAge[1],
  #                         ifelse(as.numeric(key) == 2, L0_meanAge[2],
  #                                ifelse(as.numeric(key) == 3, L0_meanAge[3],
  #                                       ifelse(as.numeric(key) == 4, L0_meanAge[4], L0_meanAge[5])))))

ggplot()+
  geom_boxplot(data = lambda_pred, aes(x=key, y=value, fill= as.factor(group_id))) +
  #scale_fill_manual(values = rep('#4890ab', 66)) + 
  guides(fill = FALSE) +  myTheme +
  scale_y_continuous(trans="log10", limits = c(3e-2, 0.3), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(title = bquote('Variation in ' ~ N(0) ~  ' within mice'),
       x = 'Mouse ID', y= bquote( ~ N(0) ~ 'count'))

ggplot()+
  geom_violin(data = lambda_pred, aes(x=key, y=value, fill= as.factor(group_id))) +
  guides(fill = FALSE) + myTheme +
  scale_y_continuous(trans="log10", limits = c(3e-2, 0.3), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(title = bquote('Variation in ' ~ lambda ~  ' within cohorts'),
       x = 'Cohort ID', y= bquote( ~ lambda(a == 0)))


ggplot()+
  geom_boxplot(data = lambda_pred, aes(x=group_id, y=value, fill= as.factor(key))) +
  scale_x_continuous(limits = c(0.7, 240), trans="log10", breaks = c(1, 3, 10, 30, 100, 300)) +
  scale_fill_discrete(name = 'Cohort') +
  scale_y_log10()+
  myTheme +
  labs(title = bquote('Variation in ' ~ lambda ~  ' with host age'),
       x = 'Host age', y= bquote( lambda(a==0)))



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


