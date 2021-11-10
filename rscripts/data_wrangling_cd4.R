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

### loading data for SP4
donor_SP4_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP4'))%>% na.omit() 


donor_SP4_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 7) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP4')) %>% na.omit() 


donor_SP4_data <- left_join(donor_SP4_ki, donor_SP4_counts, 
                            by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                            suffix = c('_ki', '_counts')) %>% na.omit() 


host_SP4_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP4'))


host_SP4_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 8) %>%
  select(contains('mouse'), contains('age'), contains('TH.SP4')) 


host_SP4_data <- left_join(host_SP4_ki, host_SP4_counts, 
                           by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                           suffix = c('_ki', '_counts')) %>% na.omit() 


SP4_data <- full_join(donor_SP4_data, host_SP4_data, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                      suffix = c('_donor', '_host'))%>%
  select(contains('mouse'), contains('age'), contains('counts'), contains('ki')) %>%
  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
  mutate(total_counts = TH.SP4_counts_donor + TH.SP4_counts_host,
         fd = TH.SP4_counts_donor/total_counts,
         ki67_donor = TH.SP4_ki_donor/100,
         ki67_host = TH.SP4_ki_host/100,
         total_ki67 = (ki67_donor * TH.SP4_counts_donor + ki67_host * TH.SP4_counts_host)/total_counts) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))%>%
  select('mouse.ID', contains('age'), 'total_counts', "fd", contains('ki67')) 

SP4_ki <- full_join(donor_SP4_data, host_SP4_data, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                    suffix = c('_donor', '_host'))%>%
  select(contains('mouse'), contains('age'), contains('ki')) %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'ki_prop')%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))

ggplot(SP4_data) +
  geom_point(aes(x= age.at.S1K, y = total_ki67), size=2) +
  geom_point(aes(x= age.at.S1K, y = ki67_donor), col=2)+
  geom_point(aes(x= age.at.S1K, y = ki67_host), col=4) + ylim(0,1)

#fitiing total dp1 numbers
sp4counts.fit <- aov(log(total_counts) ~ age.at.S1K , data = SP4_data)
summary(sp4counts.fit)
coef(sp4counts.fit)

eps_func <- function(t, a, b, c) {exp(-a * (t+b)) + c}

sp4_ki.nlm <- nls(total_ki67 ~ eps_func(age.at.S1K, a, b, c), data = SP4_data,
                  start = list(a=0.034, b =3.4, c=0.4))
summary(sp4_ki.nlm)
coef(sp4_ki.nlm)

ggplot(SP4_data)+
  geom_point(aes(x= age.at.S1K, y = total_counts))+
  geom_line(aes(x= age.at.S1K, y=exp(fitted(sp4counts.fit))), col=2, size=1.5)

ggplot(SP4_data)+
  geom_point(aes(x= age.at.S1K, y = total_ki67))+
  geom_line(aes(x= age.at.S1K, y=fitted(sp4_ki.nlm)), col=4, size=1.5)+
  ylim(0,1)

## thymic precursors 

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

write.csv(DP1_data, file = "datafiles/DP1_data.csv", row.names = F)


DP1_ki <- full_join(donor_DP1_data, host_DP1_data, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                    suffix = c('_donor', '_host'))%>%
  select(contains('mouse'), contains('age'), contains('ki')) %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'ki_prop')%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))


ggplot() +
  geom_point(data=SP4_data,aes(x = age.at.S1K, y = total_counts), size = 2) +
  geom_point(data=DP1_data,aes(x = age.at.S1K, y = total_counts), col = 2,size = 1.5) +
  #scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic SP4 cells") +
  scale_y_continuous(limits = c(1e5,3e8), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific) + myTheme +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) 

#ggsave(filename = "plots/sp4_counts.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot() +
  geom_point(data = SP4_data, aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1.5) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic SP4 cells") +
  scale_y_continuous(limits = c(1e4, 1e7), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific) + myTheme +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) 


ggplot(DP1_data) +
  geom_point(aes(x = age.at.S1K - age.at.BMT, y = fd), size = 1.5) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Days post BMT', y = NULL, title = "Donor fraction in thymic DP1 cells") + myTheme +
  scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) 

#ggsave(filename = "plots/sp4_fd.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot(SP4_data) +
  geom_point(aes(x = age.at.S1K - age.at.BMT, y = ki67_ratio, col = ageBMT_bin), size = 1.5) +
  geom_hline(yintercept = 1.0, col = 'gray50', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Days post BMT', y = NULL, title = "Donor to host ratio of Ki67+ cells in SP4 subset") +
  scale_y_log10(limits = c(0.3, 5)) +
  scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/sp4_kiratio.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot(SP4_ki) +
  geom_point(aes(x = age.at.S1K, y = ki_prop,  col = subpop), size = 1.5) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in SP4 subset") +
  ylim(0.0, 100)+
  #scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/sp4_ki67.jpg", plot = last_plot(), device = 'jpeg',
 #      width = 23, height = 8, units = 'cm')

### onotgeny data
### loading data for SP4
ont_SP4_counts <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('days'), contains('SP4')) %>%
  mutate(total_sp4 = immSP4 + matSP4) %>%
  na.omit()

### loading data for SP4
ont_SP4_ki <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('days'), contains('SP4')) %>%
  na.omit()

ont_sp4_data <- full_join(ont_SP4_counts, ont_SP4_ki, by = c('mouse.id', "age.at.S1K.days"),
                          suffix = c("_counts", "_ki")) %>%
  mutate(imm_ki = immSP4_counts * immSP4_ki /100,
         mat_ki = matSP4_counts * matSP4_ki /100,
         total_ki = (imm_ki + mat_ki)/ total_sp4)%>%
  select(contains('mouse'), contains('days'), contains('total')) 

eps_spline <- function(Time){
  eps_0 = 0.14965320; eps_f = 0.03470231; A = 3.43078629;
  exp(- eps_f * (Time + A)) + eps_0;
}

theta_spline <- function(Time){
  t0 = 1
  theta0  =  4.3E5;  theta_f = 1.8E3;  n = 2.1;   X = 30.0;   q = 3.7;
  
  theta0 + (theta_f * (Time - t0)^n) * (1 - (((Time - t0)^q)/((X^q) + ((Time - t0)^q))))
}

tseq <- seq(5, 450)
theta_vec <- sapply(tseq, theta_spline)
eps_vec <- sapply(tseq, eps_spline)


ggplot() +
  geom_point(data = SP4_data, aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1.25) +
  geom_point(data = ont_sp4_data, aes(x = age.at.S1K.days, y = total_sp4), col = 'navy', size = 2) +
 # geom_point(data = ont_sp4_data, aes(x = age.at.S1K.days, y = total_sp4), col = 4, size = 1) +
  #geom_line(aes(x = tseq, y = theta_vec), col = 'navy', size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic SP4 cells") +
  scale_y_continuous(limits = c(1e5, 1e7), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific)  +
 # scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) + 
  myTheme

#ggsave(filename = "plots/sp4_counts_overlap.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot() +
  geom_point(data = SP4_ki, aes(x = age.at.S1K, y = ki_prop,  col = subpop), size = 1.25) +
  geom_point(data = ont_sp4_data, aes(x = age.at.S1K.days, y = total_ki * 100), col = 'navy', size = 2) +
  #geom_point(data = ont_SP4_ki, aes(x = age.at.S1K.days, y = matSP4), col = 4, size = 1) +
  #geom_line(aes(x = tseq, y = eps_vec * 100), col = 'navy', size = 1) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in SP4 subset") +
  ylim(0.0, 100)+
  #scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  #scale_x_log10(limits = c(5, 450), breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/sp4_ki67_overlap.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')



### writing datafiles
#write.csv(SP4_data, file = "datafiles/sp4_data.csv", row.names = F)

#########################################################################################################
## peripheral naive T cells 

### loading data for ncd4
donor_cd4_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('age'), contains('4nai'))


donor_cd4_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 7) %>%
  select(contains('mouse'), contains('age'), contains('4nai')) 


donor_cd4_data <- left_join(donor_cd4_ki, donor_cd4_counts,
                            by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                            suffix = c('_ki', '_counts')) %>% na.omit() %>%
  mutate(spleen_ki = (SP.4nai_ki/100) * SP.4nai_counts,
         ln_ki = (LN.4nai_ki/100) * LN.4nai_counts,
         total_ki = spleen_ki + ln_ki,
         donor_counts = SP.4nai_counts + LN.4nai_counts,
         donor_ki = total_ki/donor_counts) %>%
  select(contains('mouse'), contains('age'), contains('donor'), 'total_ki') 

host_cd4_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('age'), contains('4nai'))


host_cd4_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 8) %>%
  select(contains('mouse'), contains('age'), contains('4nai')) 


host_cd4_data <- left_join(host_cd4_ki, host_cd4_counts,
                           by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                           suffix = c('_ki', '_counts')) %>% na.omit() %>%
  mutate(spleen_ki = (SP.4nai_ki/100) * SP.4nai_counts,
         ln_ki = (LN.4nai_ki/100) * LN.4nai_counts,
         total_ki = spleen_ki + ln_ki,
         host_counts = SP.4nai_counts + LN.4nai_counts,
         host_ki = total_ki/host_counts) %>%
  select(contains('mouse'), contains('age'), contains('host'), 'total_ki') 


SP4_data_fd <- SP4_data %>%
  filter(mouse.ID %in% donor_cd4_data$mouse.ID)

DP1_data_fd <- DP1_data %>%
  filter(mouse.ID %in% donor_cd4_data$mouse.ID)

#NCD4_data <- full_join(donor_cd4_data, host_cd4_data,
#                       by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
#                      suffix = c('_donor', '_host'))%>% na.omit() %>%
#  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
#  mutate(total_counts = donor_counts + host_counts,
#         total_kiprop = (total_ki_donor + total_ki_host)/total_counts,
#         fraction_donor = donor_counts/total_counts,
#         ki_dh_ratio = donor_ki/host_ki,
#         Nfd = fraction_donor/SP4_data_fd$fd) %>%
#  select(contains('age'), 'total_counts',  'donor_ki', 'host_ki', 'Nfd') %>%
#  #select(-contains('donor'), -contains('host')) %>%
#  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
#                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
#                                    'ageBMT_group3')))%>%
#  filter(Nfd <= 1.2)

NCD4_data_DP1 <- full_join(donor_cd4_data, host_cd4_data,
                       by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                       suffix = c('_donor', '_host'))%>% na.omit() %>%
  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
  mutate(total_counts = donor_counts + host_counts,
         total_kiprop = (total_ki_donor + total_ki_host)/total_counts,
         fraction_donor = donor_counts/total_counts,
         ki_dh_ratio = donor_ki/host_ki,
         Nfd = fraction_donor/DP1_data_fd$fd) %>%
  select("mouse.ID",contains('age'), 'total_counts',  'total_kiprop', 'Nfd') %>%
  #select(-contains('donor'), -contains('host')) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3'))) %>%
  filter(Nfd <= 1.2)


write.csv(NCD4_data_DP1, file = "datafiles/cd4_new.csv", row.names = F)


NCD4_ki <- full_join(donor_cd4_data, host_cd4_data,
                     by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                    suffix = c('_donor', '_host'))%>% na.omit() %>%
  select(contains('mouse'), contains('age'), 'donor_ki', 'host_ki') %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'ki_prop')%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3'))) %>% na.omit()

### writing datafiles
write.csv(NCD4_ki, file = "datafiles/cd4_ki.csv", row.names = F)



ggplot(NCD4_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of naive CD4 T cells") +
  scale_y_log10(limits = c(1e6, 1e8)) + myTheme +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) 

ggplot() +
  geom_point(data =NCD4_data2, aes(x = age.at.S1K, y = Nfd), col=4, size = 2) +
  #geom_point(data =NCD4_data2, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin),  size = 2) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Donor fraction in naive CD4 cells normalised to chimerism in DP1") +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(limits = c(60, 450), breaks = c(75, 150, 300, 450))  +
  myTheme

ggplot() +
  geom_point(data =NCD4_data, aes(x = age.at.S1K, y = Nfd), col=2, size = 2) +
  #geom_point(data =NCD4_data, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin),  size = 2) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Donor fraction in naive CD4 cells normalised to chimerism in SP4 ") +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(limits = c(60, 450), breaks = c(75, 150, 300, 450))  +
  myTheme


ggsave(filename = "plots/ncd4_Nfd.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')

ggplot(NCD4_data) +
  geom_point(data = SP4_data, aes(x = age.at.S1K - age.at.BMT, y = ki67_ratio), 
             stroke = 0.0, alpha = 0.5, shape = 19, col = 'gray40', size = 2) +
  geom_point(aes(x = age.at.S1K- age.at.BMT, y = ki_dh_ratio, col = ageBMT_bin), size = 1) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Time post BMT', y = NULL, title = "Donor to host ratio of Ki67+ cells in naive CD4 subset") +
  scale_y_log10(limits = c(0.1, 20), breaks = c(0.1, 0.3, 1, 3, 10)) +
  scale_x_log10(limits = c(10, 450), breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme 

ggsave(filename = "plots/ncd4_kiratio.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')


ggplot(NCD4_ki) +
  geom_point(aes(x = age.at.S1K, y = ki_prop * 100,  col = subpop), size = 1) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in naive CD4 subset") +
  #ylim(0.0, 100)+
  scale_y_log10(limits = c(0.1, 100)) +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme 


ggsave(filename = "plots/ncd4_kiprop.jpg", plot = last_plot(), device = 'jpeg',
       width = 23, height = 8, units = 'cm')



### loading ontogeny data for cd4s
ont_cd4_data <- read.csv("datafiles/original_data/cd4_ln.csv") 

ggplot() +
  geom_point(data = NCD4_data, aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 1) +
  geom_point(data = ont_cd4_data, aes(x = time, y = counts), col = 'navy', size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of naive CD4 cells") +
  scale_y_continuous(limits = c(1e5, 5e7), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific)  +
  # scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) + 
  myTheme

ggsave(filename = "plots/ncd4_counts_overlap.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')


ggplot() +
  #geom_point(data = NCD4_ki, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size = 2) +
  geom_point(data = NCD4_data, aes(x = age.at.S1K, y = total_kiprop * 100), size = 2) +
  geom_point(data = ont_cd4_data, aes(x = time, y = ki67 * 100), col = 6, size = 2) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in naive CD4 compartment") +
  #ylim(0.0, 100)+
  scale_y_log10(limits = c(0.3, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme 

ggsave(filename = "plots/ncd4_ki67_overlap.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')


