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

Population <- 'cd8'
ModelName <- "asm_deltavar"
#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv") 

#### Model simulations
rstan::expose_stan_functions("stan_models/CD4CRE_asm_deltavar_cd8.stan")

#Importing fitted parameters
ParamsFile <- read.csv(file.path("CD8_Chifit_output", ModelName, paste0("params_", Population, "_", ModelName, ".csv")))
params <- ParamsFile$mean[1:4]



N0vec <- c(naive_df %>% filter(group_id==1) %>%
              filter(age_anim == 14) %>% summarise(mean(naive_labelled)),
            naive_df %>% filter(group_id==2) %>%
              filter(age_anim == 28) %>% summarise(mean(naive_labelled)),
            naive_df %>% filter(group_id==3) %>%
              filter(age_anim == 42) %>% summarise(mean(naive_labelled)),
            naive_df %>% filter(group_id==4) %>%
              filter(age_anim == 70) %>% summarise(mean(naive_labelled)),
            naive_df %>% filter(group_id==5) %>%
              filter(age_anim == 189) %>% summarise(mean(naive_labelled)))

N0_vec <- unlist(N0vec)

## prediction points
tseq1 <- seq(14, 120, length.out = 100)
tseq2 <- seq(28, 120, length.out = 100)
tseq3 <- seq(42, 100, length.out = 100)
tseq4 <- seq(70, 160, length.out = 100)
tseq5 <- seq(189, 280, length.out = 100)


parstan1_vec <- c(N0_vec[1], params[2:4], c(1, 14))
parstan2_vec <- c(N0_vec[2], params[2:4], c(7, 28))
parstan3_vec <- c(N0_vec[3], params[2:4], c(28, 42))
parstan4_vec <- c(N0_vec[4], params[2:4], c(56, 70))
parstan5_vec <- c(N0_vec[5], params[2:4], c(175, 189))


## stan solution
sim1 <- sapply(tseq1,  N_total_time, parms = parstan1_vec)
sim2 <- sapply(tseq2,  N_total_time, parms = parstan2_vec)
sim3 <- sapply(tseq3,  N_total_time, parms = parstan3_vec)
sim4 <- sapply(tseq4,  N_total_time, parms = parstan4_vec)
sim5 <- sapply(tseq5,  N_total_time, parms = parstan5_vec)

## labels for individual age cohorts
age_cohorts <- c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100), rep(5, 100))
## predictions data frame
naive_sim_df <- data.frame("timeseries" = c(tseq1, tseq2, tseq3, tseq4, tseq5),
                           "group_id" = age_cohorts,
                           "median" = c(sim1, sim2, sim3, sim4, sim5))


ggplot()+
  geom_line(data = naive_sim_df, aes(x = timeseries, y = median, col = as.factor(group_id)),
             size =1.2) +
  geom_hline(yintercept = N0_vec[5], size=2, col='darkred')+
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(group_id)), size =3) +
  #geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled), col = "black", size =1.2) +
  scale_color_discrete(name = "Group ID", guide = guide_legend(nrow = 7)) +
  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +guides(col = F) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
  myTheme 











