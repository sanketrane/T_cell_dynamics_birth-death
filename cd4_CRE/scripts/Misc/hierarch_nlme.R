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
  select('id', 'age_anim', 'age_group', 'naive_labelled') %>%
  filter(age_group == 1) %>%
  rename(mouse_id = id) %>%
  mutate(host_age = age_anim -14)


ggplot()+
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                 col = as.factor(mouse_id)), size =3) +
  geom_line(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(mouse_id)), size =1.2) +
  scale_color_discrete(name = "Mouse ID", guide = guide_legend(nrow = 3)) +
  scale_y_log10(limits = c(5e3, 1e6), breaks = c(1e5, 1e6, 1e4, 1e7)) +
  scale_x_continuous(limits = c(0, 120), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_rect(fill = NA, colour = NA))



library(nlme)
datasetG <- groupedData(naive_labelled ~ age_anim | mouse_id, naive_df)

f1 <- log(naive_labelled) ~ 10^N0 * exp(- lambda * age_anim)
nlin.mix <- nlme(f1,
                data = datasetG,
                fixed = list(N0 ~ 1, lambda ~ 1),
                random = lambda ~ 1 | mouse_id,
                start= c(N0 = 5.6, lambda = 0.05))

summary(nlin.mix)
plot(nlin.mix)

