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

### loading data for DP1
donor_DP1_counts <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('age'), contains('TH.DP1'))


donor_DP1_ki <- read_excel("datafiles/UPDATED_master_doc.xlsx", sheet = 7) %>%
  select(contains('mouse'), contains('age'), contains('TH.DP1')) 


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


DP1_data <- full_join(donor_DP1_data, host_DP1_data, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                      suffix = c('_donor', '_host'))%>%
  select(contains('mouse'), contains('age'), contains('counts'), contains('ki')) %>%
  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
  mutate(total_counts = TH.DP1_counts_donor + TH.DP1_counts_host,
         fd = TH.DP1_counts_donor/total_counts,
         ki67_donor = TH.DP1_ki_donor/100,
         ki67_host = TH.DP1_ki_host/100,
         ki67_ratio = ki67_donor/ki67_host) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))%>%
  select('mouse.ID', contains('age'), 'total_counts', "fd", contains('ki67')) 

DP1_ki <- full_join(donor_DP1_data, host_DP1_data, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                    suffix = c('_donor', '_host'))%>%
  select(contains('mouse'), contains('age'), contains('ki')) %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'ki_prop')%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))


ggplot(DP1_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 2) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic DP1 cells") +
  scale_y_continuous(limits = c(1e6, 2e8), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific) + myTheme +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) 

#ggsave(filename = "plots/DP1_counts.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')



ggplot(DP1_data) +
  geom_point(aes(x = age.at.S1K - age.at.BMT, y = fd, col = ageBMT_bin), size = 1.5) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Days post BMT', y = NULL, title = "Donor fraction in thymic DP1 cells") + myTheme +
  scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) 

#ggsave(filename = "plots/DP1_fd.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot(DP1_data) +
  geom_point(aes(x = age.at.S1K - age.at.BMT, y = ki67_ratio, col = ageBMT_bin), size = 1.5) +
  geom_hline(yintercept = 1.0, col = 'gray50', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Days post BMT', y = NULL, title = "Donor to host ratio of Ki67+ cells in DP1 subset") +
  scale_y_log10(limits = c(0.3, 5)) +
  scale_x_log10(limits = c(10, 450), breaks = c(10, 30, 100, 300)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/DP1_kiratio.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot(DP1_ki) +
  geom_point(aes(x = age.at.S1K, y = ki_prop,  col = subpop), size = 1.5) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in DP1 subset") +
  ylim(0.0, 100)+
  #scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/DP1_ki67.jpg", plot = last_plot(), device = 'jpeg',
 #      width = 23, height = 8, units = 'cm')

### onotgeny data
### loading data for DP1
ont_DP1_counts <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 2) %>%
  select(contains('mouse'), contains('days'), contains('DP1')) %>%
  #mutate(total_DP1 = immDP1 + matDP1) %>%
  na.omit()

### loading data for DP1
ont_DP1_ki <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 3) %>%
  select(contains('mouse'), contains('days'), contains('DP1')) %>%
  na.omit()

ont_DP1_data <- full_join(ont_DP1_counts, ont_DP1_ki, by = c('mouse.id', "age.at.S1K.days"),
                          suffix = c("_counts", "_ki")) %>%
  mutate(total_ki = (DP1_ki)/DP1_counts)%>%
  select(contains('mouse'), contains('days'), contains('DP1')) 

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
  geom_point(data = DP1_data, aes(x = age.at.S1K, y = total_counts, col = ageBMT_bin), size = 2) +
  geom_point(data = ont_DP1_data, aes(x = age.at.S1K.days, y = DP1_counts), col = 'navy', size = 2) +
 # geom_point(data = ont_DP1_data, aes(x = age.at.S1K.days, y = total_DP1), col = 4, size = 1) +
  #geom_line(aes(x = tseq, y = theta_vec), col = 'navy', size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Total counts of thymic DP1 cells") +
  scale_y_continuous(limits = c(1e6, 2e8), trans = "log10",
                     minor_breaks = log10minorbreaks, labels =fancy_scientific)  +
 # scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) + 
  myTheme

#ggsave(filename = "plots/DP1_counts_overlap.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')


ggplot() +
  geom_point(data = DP1_ki, aes(x = age.at.S1K, y = ki_prop,  col = subpop), size = 2) +
  geom_point(data = ont_DP1_data, aes(x = age.at.S1K.days, y = DP1_ki), col = 'navy', size = 2) +
  #geom_point(data = ont_DP1_ki, aes(x = age.at.S1K.days, y = matDP1), col = 4, size = 1) +
  #geom_line(aes(x = tseq, y = eps_vec * 100), col = 'navy', size = 1) +
  scale_color_discrete(name = "Subset", labels = c("Donor", "Host")) +
  labs(x = 'Host age', y = NULL, title = "% Ki67high cells in DP1 subset") +
  ylim(0.0, 100)+
  #scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 0.3, 1, 3, 10, 30, 100)) +
  #scale_x_log10(limits = c(5, 450), breaks = c(75, 150, 300, 450)) +
  #facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme

#ggsave(filename = "plots/DP1_ki67_overlap.jpg", plot = last_plot(), device = 'jpeg',
#       width = 13, height = 8, units = 'cm')



### writing datafiles
#write.csv(DP1_data, file = "datafiles/DP1_data.csv", row.names = F)

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


NCD4_data <- full_join(donor_cd4_data, host_cd4_data,
                       by = c("mouse.ID", "age.at.S1K", "age.at.BMT"),
                      suffix = c('_donor', '_host'))%>% na.omit() %>%
  #gather(-c(mouse.ID, age.at.S1K, age.at.BMT), key = 'subpop', value = 'counts') %>%
  mutate(total_counts = donor_counts + host_counts,
         total_kiprop = (total_ki_donor + total_ki_host)/total_counts,
         fd = donor_counts/total_counts,
         ki_dh_ratio = donor_ki/host_ki) %>%
  select("mouse.ID", contains('age'), 'total_counts',  'donor_ki', 'host_ki', "fd") %>%
  #select(-contains('donor'), -contains('host')) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                             ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                    'ageBMT_group3')))

NCD4_fd <- NCD4_data %>%
  filter(mouse.ID %in% DP1_data$mouse.ID)%>%
  mutate(Nfd= fd/DP1_data$fd)

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

ggplot(NCD4_data) +
  geom_point(aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size = 1) +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size = 1) +
  scale_color_discrete(name = "Age at BMT", labels = c("7-9wks", "9-11wks", "11-25wks")) +
  labs(x = 'Host age', y = NULL, title = "Normalised chimerism in naive CD4 cells") + myTheme +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_log10(limits = c(60, 450), breaks = c(75, 150, 300, 450)) 


ggsave(filename = "plots/ncd4_Nfd.jpg", plot = last_plot(), device = 'jpeg',
       width = 13, height = 8, units = 'cm')

ggplot(NCD4_data) +
  geom_point(data = DP1_data, aes(x = age.at.S1K - age.at.BMT, y = ki67_ratio), 
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
ont_cd4_data <- read.csv("datafiles/cd4_ln.csv") 

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


