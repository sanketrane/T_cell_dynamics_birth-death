setwd("~/Desktop/GIt_repos/magnum_opus")

library(readxl)
library(tidyverse)

theme_set(theme_bw())
myTheme <- theme(text = element_text(size = 11), axis.text = element_text(size = 11), axis.title =  element_text(size = 10, face = "bold"),
                 plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())


rfp_bg <- read_excel("data/timestamp/RFP_bg.xlsx", sheet = 1) %>%
  gather(key = age_cell, value = rfp_bg) %>% na.omit() %>%
  mutate(host_age = ifelse(age_cell == "d14", 14,
                           ifelse(age_cell == "d28", 28,
                                  ifelse(age_cell == "d56", 56,
                                         ifelse(age_cell == "d84", 84,245)))))
#group_by(age_cell) #%>%
#summarise(mean_counts = mean(counts))


total_counts <- read_excel("data/timestamp/total_counts.xlsx", sheet = 1) %>%
  gather(key = age_cell, value = counts) %>% na.omit() %>%
  mutate(host_age = ifelse(age_cell == "d14", 14,
                           ifelse(age_cell == "d28", 28,
                                  ifelse(age_cell == "d56", 56,
                                         ifelse(age_cell == "d84", 84,245)))))
  #group_by(age_cell) #%>%
  #summarise(mean_counts = mean(counts))


rfp_func <- function(t, p, b0, b){
  
  p * b0 * exp(b*t)/(p + b0*(exp(b*t) - 1))
}

rfp_fit <- nls(rfp_bg ~ rfp_func(host_age, p, b0, b), data = rfp_bg,
               start = list(p = 0.2, b0 = 0.01, b = 0.01))

summary(rfp_fit)

counts_func <- function(t, k, n0, g){
  
    10^(k) * 10^(n0) * exp(g*t)/(10^(k) + 10^(n0) * (exp(g*t) - 1))
}

counts_func1 <- function(t, b0, r1, r2, r3){
  t0 = 14
  t1 =56
  t2 = 84
  b_t1 = 10^b0 * exp(10^r1 * (t1 - t0))
  b_t2 = b_t1 * exp(-10^r2 * (t2 - t1))
  
  value = ifelse(t < t1,
                 10^b0 * exp(10^r1 * (t - t0)),
                 ifelse(t <= t2, b_t1 * exp(-10^r2 * (t-t1)),
                        b_t2 * exp(-10^r3 * (t-t2))))
  
  return(value)
}

counts_func2 <- function(t, b0, r1, r2){
  t1 =56
  b_t1 <- 10^b0 * (1 - exp(-10^r1 * t1))
  
  value = ifelse(t < t1,
                 10^b0 * (1 - exp(-10^r1 * t)),
                 b_t1 * exp(-10^r2 * (t-t1)))
  
  return(value)
}


counts_func1(14, counts_fit2$par[1], counts_fit2$par[2],
             counts_fit2$par[3], counts_fit2$par[4])

ts_p <-  10^seq(0, 2.4, length.out = 300)
counts_vec <- counts_func1(ts_p, counts_fit2$par[1], counts_fit2$par[2],
                           counts_fit2$par[3], counts_fit2$par[4])

ggplot()+
  geom_point(data = total_counts, aes(x = host_age, y = counts), size =2) +
  geom_line(aes(x= ts_p, y = counts_vec), size =1.25) + scale_y_log10()


counts_vec <- counts_func1(ts_p, counts_fit2$par[1], counts_fit2$par[2],
                 counts_fit2$par[3], counts_fit2$par[4])

optim_func <- function(pars, boot_data){
  b0 = pars[1]  
  r1 = pars[2]  
  r2 = pars[3]  
  r3 = pars[4]  
  
  ssr = sum((log(boot_data$counts) - log(counts_func1(boot_data$host_age, b0, r1, r2, r3)))^2)
  
  k = length(pars)
  n1 = nrow(boot_data)
  
  logL = - (n1/2) * log(ssr)
  
  return(-logL)
}

counts_fit2 <- optim(par = c(b0 = 7.5, r1 = -1.8, r2 = -1.7, r3 = -3.7), fn = optim_func,
                     boot_data = total_counts, hessian = T, control = list(trace =T, maxit = 5000))

10^counts_fit2$par
ts_p <- seq(10, 300)


## importing ontogeny data from seddon lab
DataFile <- file.path("data", "cd8_ln.csv")
data_fit <- read.csv(DataFile) %>%
  mutate(age.at.S1K = time,
         time_post_t0 = age.at.S1K - min(age.at.S1K)) %>% arrange(age.at.S1K) %>% unique() %>%
  select("age.at.S1K", contains('counts'), contains('ki67'))


ggplot() +
  geom_line(aes(x= ts_p, y = counts_func(ts_p, 7.06, 5.54431, 0.15014)), size = 1.5) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = counts),  col = 2, size = 3) +
  geom_point(data = total_counts, aes(x = host_age, y = counts),  col = 4, bg=4, pch=24, size =3) +
  #geom_line(aes(x= ts_p, y = counts_func1(ts_p, counts_fit2$par[1], counts_fit2$par[2],
   #                                      counts_fit2$par[3], counts_fit2$par[4])), col = 3, size =1.25) +
  labs(title= paste0('Total numbers of naive CD8 T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 300), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 3e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8)) 

ggplot()+
  geom_point(data = rfp_bg, aes(x = host_age, y = rfp_bg), size =2) + 
  geom_line(aes(x= ts_p, y = rfp_func(ts_p, 0.19455, 0.01647, 0.09148)), size =1.25) #+
  #scale_y_log10() #+ scale_x_log10(breaks = c(20, 60, 200, 600))

ggplot()+
  geom_point(data = total_counts, aes(x = host_age, y = counts), size =2) +
  geom_line(aes(x= ts_p, y = counts_func(ts_p, 7.06867, 5.54431, 0.15014)), size =1.25) +
  scale_y_log10(limits = c(5e5, 5e7)) #+ scale_x_log10(breaks = c(20, 60, 200, 600))



tstmp <- read_excel("data/timestamp/data_03.xlsx", sheet = "data") %>%
  mutate(labelled_counts = (((fraction - rfp_func(age_anim, 0.19455, 0.01647, 0.09148))/100) 
                            * counts_func(age_anim, 7.06867, 5.54431, 0.15014)))

log10minorbreaks=as.numeric(1:10 %o% 10^(-2:8))

ggplot()+
  geom_point(data = tstmp,  aes(x = age_anim, y = fraction, col = as.factor(age_group)), alpha = 0.7, size =2) +
  scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100), 
               minor_breaks = log10minorbreaks, labels = fancy_scientific)+
  scale_color_discrete(name = "Age at timestamp (d)", guide = guide_legend(nrow = 3)) +
  labs(x = "Mouse age", title = "Counts of timestamped CD8 cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.88, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())

ggplot()+
  geom_point(data = tstmp,  aes(x = age_anim, y = labelled_counts, col = as.factor(age_group)), size =2) +
  scale_y_log10(limits = c(1e4, 1e7), breaks = c(1e5, 1e6, 1e4, 1e7), 
                minor_breaks = log10minorbreaks, labels = fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  scale_color_discrete(name = "Age at timestamp (d)", guide = guide_legend(nrow = 3)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped CD8 cells", y=NULL)  + myTheme +
  theme(legend.position = c(0.88, 0.84), legend.key.size = unit(0.6, "cm"),
        strip.background = element_blank())

























