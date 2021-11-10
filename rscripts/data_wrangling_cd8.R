###
rm(list = ls()); gc()

### loading the packages
library(readxl)
library(tidyverse)


### theme set for plots
theme_set(theme_bw())
myTheme <- theme(text = element_text(size = 10), axis.text = element_text(size = 10), 
                 axis.title =  element_text(size = 10, face = "bold"),
                 plot.title = element_text(size=10,  hjust = 0.5, face = "bold"),
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


ageBMT_names <- c(`ageBMT_group1` = "7-9wks",
                  `ageBMT_group2` = "9-11wks", 
                  `ageBMT_group3` = "11-25wks")


#########################################################################################################
## thymic precursors 

### loading data for SP8
donor_SP8_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP8'))


donor_SP8_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 7) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP8')) 


donor_SP8_data <- left_join(donor_SP8_ki, donor_SP8_counts, 
                            by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                            suffix = c('_ki', '_counts')) %>% na.omit() 


host_SP8_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP8'))


host_SP8_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 8) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP8')) 


host_SP8_data <- left_join(host_SP8_ki, host_SP8_counts, 
                           by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                           suffix = c('_ki', '_counts')) %>% na.omit() 


SP8_data <- full_join(donor_SP8_data, host_SP8_data, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                      suffix = c('_donor', '_host'))%>% 
  select(contains('mouse'), contains('age'), contains('counts'), contains('ki')) %>%
  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
  mutate(total_counts = TH.SP8_counts_donor + TH.SP8_counts_host,
         fd = TH.SP8_counts_donor/total_counts,
         ki67_donor = TH.SP8_ki_donor/100,
         ki67_host = TH.SP8_ki_host/100,
         total_ki67 = (ki67_donor * TH.SP8_counts_donor + ki67_host * TH.SP8_counts_host)/total_counts) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))%>%
  select('mouse.ID', contains('age'), 'total_counts', "fd", contains('ki67')) %>% na.omit()


SP8_ki <- full_join(donor_SP8_data, host_SP8_data, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                    suffix = c('_donor', '_host'))%>%
  select(contains('mouse'), contains('age'), contains('ki')) %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'ki_prop')%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3'))) %>% na.omit()


#fitiing total dp1 numbers
sp8counts.fit <- aov(log(total_counts) ~ age.at.S1K , data = SP8_data)
summary(sp8counts.fit)
coef(sp8counts.fit)

eps_func <- function(t, a, b, c) {exp(-a * (t+b)) + c}

sp8_ki.nlm <- nls(total_ki67 ~ eps_func(age.at.S1K, a, b, c), data = SP8_data,
                  start = list(a=0.034, b =3.4, c=0.4))
summary(sp8_ki.nlm)
coef(sp8_ki.nlm)

ggplot(SP8_data)+
  geom_point(aes(x= age.at.S1K, y = total_counts))+
  geom_line(aes(x= age.at.S1K, y=exp(fitted(sp8counts.fit))), col=2, size=1.5)

ggplot(SP8_data)+
  geom_point(aes(x= age.at.S1K, y = total_ki67))+
  geom_line(aes(x= age.at.S1K, y=fitted(sp8_ki.nlm)), col=4, size=1.5)+
  ylim(0,1)


### loading data for DP1
donor_DP1_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('age'), contains('TH.DP1'))%>% na.omit() 


donor_DP1_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 7) %>%
  select(contains('mouse'), contains('age'), contains('TH.DP1')) %>% na.omit() 


donor_DP1_data <- left_join(donor_DP1_ki, donor_DP1_counts, 
                            by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                            suffix = c('_ki', '_counts')) %>% na.omit() 


host_DP1_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('age'), contains('TH.DP1'))


host_DP1_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 8) %>%
  select(contains('mouse'), contains('age'), contains('TH.DP1')) 


host_DP1_data <- left_join(host_DP1_ki, host_DP1_counts, 
                           by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                           suffix = c('_ki', '_counts')) %>% na.omit() 


DP1_data <- full_join(donor_DP1_counts, host_DP1_counts, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                      suffix = c('_donor', '_host'))%>%
  select(contains('mouse'), contains('age'), contains('DP1')) %>%
  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
  mutate(total_counts = TH.DP1_donor + TH.DP1_host,
         fd = TH.DP1_donor/total_counts) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))%>%
  select('mouse.ID', contains('age'), 'total_counts', "fd", contains('ki67')) %>%
  na.omit()



ggplot(SP8_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1.5) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic SP8 cells") +
  scale_y_continuous(limits = c(1e4, 1e7), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific) + myTheme +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) 

#ggsave(filename = "plots/SP8_counts.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot(SP8_data) +
  geom_point(aes(x = age.at.S1K - age.at.BMT, y = fd, col = ageBMT_bin), size = 1.5) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Days post BMT', y = NULL, title = "Donor fraction in thymic SP8 cells") + myTheme +
  scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) 

#ggsave(filename = "plots/SP8_fd.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot(SP8_data) +
  geom_point(aes(x = age.at.S1K - age.at.BMT, y = ki67_ratio, col = ageBMT_bin), size = 1.5) +
  geom_hline(yintercept = 1.0, col = 'gray50', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Days post BMT', y = NULL, title = "Donor to host ratio of Ki67+ cells in SP8 subset") +
  scale_y_log10(limits = c(0.3, 5)) +
  scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/SP8_kiratio.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot(SP8_ki) +
  geom_point(aes(x = age.at.S1K - age.at.BMT, y = ki_prop,  col = subpop), size = 1.5) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in SP8 subset") +
  ylim(0.0, 100)+
  #scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  scale_x_log10(limits = c(10, 450), breaks = c(75, 150, 300, 450)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/SP8_ki67.jpg", plot = last_plot(), device = 'jpeg',
#       width = 23, height = 8, units = 'cm')

### onotgeny data
### loading data for SP8
ont_SP8_counts <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('days'), contains('SP8')) %>%
  mutate(total_SP8 = immSP8 + matSP8) %>%
  na.omit()

### loading data for SP4
ont_SP8_ki <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('days'), contains('SP8')) %>%
  na.omit()

ont_SP8_data <- full_join(ont_SP8_counts, ont_SP8_ki, by = c('mouse.id', "age.at.S1K.days"),
                          suffix = c("_counts", "_ki")) %>%
  mutate(imm_ki = immSP8_counts * immSP8_ki /100,
         mat_ki = matSP8_counts * matSP8_ki /100,
         total_ki = (imm_ki + mat_ki)/ total_SP8)%>%
  select(contains('mouse'), contains('days'), contains('mat')) 

eps_spline <- function(Time){
  eps_0 = 0.24510453; eps_f = 0.01559996; A = 14.83715328;
  exp(- eps_f * (Time + A)) + eps_0;
}

theta_spline2 <- function(Time){
  t0 = 1
  theta0  =  9E4;  theta_f = 68;  n = 3;   X = 25;   q = 3.75;
  
  theta0 + (theta_f * (Time - t0)^n) * (1 - (((Time - t0)^q)/((X^q) + ((Time - t0)^q))))
}

tseq <- seq(5, 450)
theta_vec <- sapply(tseq, theta_spline)
eps_vec <- sapply(tseq, eps_spline)

chivec <- sapply(tseq - 40, Chi_spline)

donortheta_vec <- sapply(tseq, theta_donor, parms=parstan_vec)
hosttheta_vec <- sapply(tseq, theta_host, parms=parstan_vec)
total_vec <- sapply(tseq, theta_spline, parms=parstan_vec)

ggplot() +
  geom_point(data = SP8_data, aes(x = age.at.S1K, y = total_counts * fd), col=2, size = 2) +
  geom_point(data = SP8_data, aes(x = age.at.S1K, y = total_counts * (1 - fd)), col=4, size = 2) +
  geom_line(aes(x = tseq, y = donortheta_vec ), col = 2, size = 1) +
  geom_line(aes(x = tseq, y = hosttheta_vec), col = 4, size = 1) +
  geom_line(aes(x = tseq, y = total_vec), col = 1, size = 1) +
  geom_line(aes(x = tseq, y = donortheta_vec + hosttheta_vec), col = 3, size = 1) +
  #geom_point(data = SP8_data, aes(x = age.at.S1K, y = total_counts), col='navy', size = 2) +
  #geom_point(data = ont_SP8_data, aes(x = age.at.S1K.days, y = matSP8_counts), col = 3, size = 2) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic SP8 cells") +
  scale_y_continuous(limits = c(1e3, 5e6), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific)  +
  scale_x_log10(limits = c(10, 450), breaks = c(75, 150, 300, 450)) + 
  myTheme


ggplot() +
  geom_point(data = SP8_data, aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1.25) +
  geom_point(data = ont_SP8_data, aes(x = age.at.S1K.days, y = matSP8_counts), col = 'navy', size = 2) +
  # geom_point(data = ont_SP8_data, aes(x = age.at.S1K.days, y = total_SP8), col = 4, size = 1) +
  geom_line(aes(x = tseq, y = theta_vec), col = 'navy', size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic SP8 cells") +
  scale_y_continuous(limits = c(1e4, 5e6), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific)  +
  # scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) + 
  myTheme

#ggsave(filename = "plots/SP8_counts_overlap.jpg", plot = last_plot(), device = 'jpeg',
 #      width = 13, height = 8, units = 'cm')


ggplot() +
  geom_point(data = SP8_ki, aes(x = age.at.S1K, y = ki_prop,  col = subpop), size = 1.25) +
  geom_point(data = ont_SP8_data, aes(x = age.at.S1K.days, y = matSP8_ki), col = 'navy', size = 2) +
  #geom_point(data = ont_SP8_ki, aes(x = age.at.S1K.days, y = matSP8), col = 4, size = 1) +
  geom_line(aes(x = tseq, y = eps_vec * 100), col = 'navy', size = 1) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in SP8 subset") +
  ylim(0.0, 100)+
  #scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  #scale_x_log10(limits = c(5, 450), breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme


### writing datafiles
write.csv(SP8_data, file = "datafiles/sp8_data.csv", row.names = F)

#########################################################################################################
## peripheral naive T cells 

### loading data for ncd8
donor_cd8_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('age'), contains('8nai'))


donor_cd8_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 7) %>%
  select(contains('mouse'), contains('age'), contains('8nai')) 


donor_cd8_data <- left_join(donor_cd8_ki, donor_cd8_counts,
                            by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                            suffix = c('_ki', '_counts')) %>% na.omit() %>%
  mutate(spleen_ki = (SP.8nai_ki/100) * SP.8nai_counts,
         ln_ki = (LN.8nai_ki/100) * LN.8nai_counts,
         total_ki = spleen_ki + ln_ki,
         donor_counts = SP.8nai_counts + LN.8nai_counts,
         donor_ki = total_ki/donor_counts) %>%
  select(contains('mouse'), contains('age'), contains('donor'), 'total_ki') 

host_cd8_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('age'), contains('8nai'))


host_cd8_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 8) %>%
  select(contains('mouse'), contains('age'), contains('8nai')) 


host_cd8_data <- left_join(host_cd8_ki, host_cd8_counts,
                           by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                           suffix = c('_ki', '_counts')) %>% na.omit() %>%
  mutate(spleen_ki = (SP.8nai_ki/100) * SP.8nai_counts,
         ln_ki = (LN.8nai_ki/100) * LN.8nai_counts,
         total_ki = spleen_ki + ln_ki,
         host_counts = SP.8nai_counts + LN.8nai_counts,
         host_ki = total_ki/host_counts) %>%
  select(contains('mouse'), contains('age'), contains('host'), 'total_ki') 


SP8_data_fd <- SP8_data %>%
  filter(mouse.ID %in% donor_cd8_data$mouse.ID)


DP1_data_fd <- DP1_data %>%
  filter(mouse.ID %in% donor_cd8_data$mouse.ID)

#NCD8_data <- full_join(donor_cd8_data, host_cd8_data,
#                       by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
#                      suffix = c('_donor', '_host')) %>%
#  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
#  mutate(total_counts = donor_counts + host_counts,
#         total_kiprop = (total_ki_donor + total_ki_host)/total_counts,
#         fd = donor_counts/total_counts,
#         ki_dh_ratio = donor_ki/host_ki,
#         Nfd = fd/SP8_data_fd$fd) %>%
#  select(contains('age'), 'total_counts',  'donor_ki', 'host_ki', 'Nfd') %>%
#  #select(-contains('donor'), -contains('host')) %>%
#  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
#                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
#                                    'ageBMT_group3'))) %>% na.omit()%>%
#  filter(Nfd <= 1.2)

NCD8_data_DP1 <- full_join(donor_cd8_data, host_cd8_data,
                       by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                       suffix = c('_donor', '_host')) %>%
  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
  mutate(total_counts = donor_counts + host_counts,
         total_kiprop = (total_ki_donor + total_ki_host)/total_counts,
         fd = donor_counts/total_counts,
         ki_dh_ratio = donor_ki/host_ki,
         Nfd = fd/DP1_data_fd$fd) %>%
  select("mouse.ID", contains('age'), 'total_counts',  'total_kiprop', 'Nfd') %>%
  #select(-contains('donor'), -contains('host')) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3'))) %>% na.omit()%>%
  filter(Nfd <= 1.2)


### writing datafiles
write.csv(NCD8_data_DP1, file = "datafiles/cd8_new.csv", row.names = F)

NCD8_ki <- full_join(donor_cd8_data, host_cd8_data,
                     by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                    suffix = c('_donor', '_host')) %>%
  select(contains('mouse'), contains('age'), 'donor_ki', 'host_ki') %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'ki_prop')%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3'))) %>% na.omit()

### writing datafiles
write.csv(NCD8_ki, file = "datafiles/cd8_ki.csv", row.names = F)


ggplot(NCD8_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of naive CD8 T cells") +
  scale_y_continuous(limits = c(1e6, 1e8), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific) + myTheme +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) 

ggsave(filename = "plots/ncd8_counts.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')


ggplot(NCD8_data) +
  #geom_point(aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size = 1) +
  geom_point(aes(x = age.at.S1K, y = Nfd), col=2, size = 1) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Normalised chimerism in naive CD8 cells") + 
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(limits = c(60, 450), trans = 'log10', breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme 

ggplot(NCD8_data2) +
  #geom_point(aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size = 1) +
  geom_point(aes(x = age.at.S1K, y = Nfd), col=4,  size = 1) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Normalised chimerism in naive CD8 cells") + 
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(limits = c(60, 450), trans = 'log10', breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme 

ggsave(filename = "plots/ncd8_Nfd.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')


ggplot(NCD8_data) +
  geom_point(data = SP8_data, aes(x = age.at.S1K - age.at.BMT, y = ki67_ratio), 
             stroke = 0.0, alpha = 0.5, shape = 19, col = 'gray40', size = 2) +
  geom_point(aes(x = age.at.S1K- age.at.BMT, y = ki_dh_ratio, col = ageBMT_bin), size = 1) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Time post BMT', y = NULL, title = "Donor to host ratio of Ki67+ cells in naive CD8 subset") +
  scale_y_log10(limits = c(0.1, 20), breaks = c(0.1, 0.3, 1, 3, 10)) +
  scale_x_log10(limits = c(10, 450), breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme 

ggsave(filename = "plots/ncd8_kiratio.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')


ggplot(NCD8_ki) +
  geom_point(aes(x = age.at.S1K, y = ki_prop * 100,  col = subpop), size = 1) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67+ cells in naive CD8 subset") +
  scale_y_log10(limits = c(0.1, 100)) + myTheme +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names))


ggsave(filename = "plots/ncd8_kiprop.jpg", plot = last_plot(), device = 'jpeg',
       width = 23, height = 8, units = 'cm')


### loading ontogeny data for cd8s
ont_cd8_data <- read.csv("datafiles/cd8_ln.csv") 

ggplot() +
  geom_point(data = NCD8_data, aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1) +
  geom_point(data = ont_cd8_data, aes(x = time, y = counts), col = 'navy', size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of naive cd8 cells") +
  scale_y_continuous(limits = c(1e5, 5e7), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific)  +
  # scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) + 
  myTheme

ggsave(filename = "plots/ncd8_counts_overlap.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')


ggplot() +
  geom_point(data = NCD8_ki, aes(x = age.at.S1K, y = ki_prop * 100,  col = subpop), size = 1) +
  geom_point(data = ont_cd8_data, aes(x = time, y = ki67 * 100), col = 'navy', size = 1) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in naive cd8 compartment") +
  #ylim(0.0, 100)+
  scale_y_log10(limits = c(0.3, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  #scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme 

ggsave(filename = "plots/ncd8_ki67_overlap.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')





