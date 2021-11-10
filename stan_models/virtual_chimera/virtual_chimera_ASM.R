rm(list = ls()); gc()
library(rstan)
library(tidyverse)

theme_set(theme_bw())
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), 
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
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

substrRight <- function(x, n){
  substr(x, 1, nchar(x)-n)
}


ageBMT_names <- c(`ageBMT_group1` = "7-9wks",
                  `ageBMT_group2` = "9-11wks", 
                  `ageBMT_group3` = "11-25wks")


## importing ASM stan model
ModelName <- 'asmS52_rhovar'
Population <- "cd8"
OutputDir <- file.path("output_csv")
SaveDir <- file.path('out_fit', Population, ModelName) 

stan_file <- paste0("asm_rhovar", '_', Population, '.stan')
expose_stan_functions(file.path('stan_models/virtual_chimera', stan_file))


## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) %>%
  arrange(age.at.S1K)

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34))%>%
  arrange(age.at.S1K)

data_counts <- chimera_data %>% 
  select(contains('age'), contains('counts'), contains('kiprop')) %>%
  mutate(dataset = ifelse(age.at.BMT <= 63, 'ageBMT_group1',
                          ifelse(age.at.BMT <= 77, 'ageBMT_group2',
                                 'ageBMT_group3'))) %>%
  bind_rows(ontogeny_data) %>%
  arrange(age.at.S1K)

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file) 



ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_counts, col = dataset), size =2.5) + 
  labs(x = 'Host age', y = NULL, title = paste0("Total counts of naive ", Population, " T cells")) +
  scale_y_log10(limits = c(1e5, 1e8)) 

ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_kiprop * 100, col = dataset), size =2.5) + 
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67+ cells in naive ", Population, " T cells")) +
  scale_y_log10(limits = c(0.3, 100)) 

ggplot() +
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = ki_dh_ratio, col = ageBMT_bin), size =2.5) + 
  labs(x = 'Host age', y = NULL, title = paste0("Ratio of donor to host Ki67+ cell in naive ", Population, " T cells")) +
  scale_y_log10(limits = c(0.3, 100)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + guides(col = F)



# vector of parameters
param_file <- file.path(SaveDir, paste0('params_', Population, '_', ModelName, '.csv'))
par_df <- data.frame(read.csv(param_file))[1:6, 'mean']
parstan_vec <- c(par_df[1:4]) 

## time course for predictions
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)
tb_pred1 <- rep(54, 300)
tb_pred2 <- rep(71, 300)
tb_pred3 <- rep(97, 300)

data_ont <- data_counts$age.at.S1K                        # timepoints in observations 
data_chi <- chimera_data$age.at.S1K 
tb_chi <- chimera_data$age.at.BMT

ntotal_counts <- N_total_time(data_ont, parstan_vec)
npool_counts <- N_total_time(data_chi, parstan_vec)
nhost_counts <- N_host_time(data_chi, tb_chi, parstan_vec)
ndonor_counts <- N_donor_time(data_chi, tb_chi, parstan_vec)
chivec <- sapply(data_chi, Chi_spline)


v1 <- Sys.time()
kpos_vec <- U_total_time(data_ont, parstan_vec)
kpos_pred <- kpos_vec/ntotal_counts
v2 <- Sys.time()
v2-v1
kpos_host_counts <- U_host_time(data_chi, tb_chi, parstan_vec)
kpos_host_ki <- kpos_host_counts/nhost_counts
kpos_donor_counts <- U_donor_time(data_chi, tb_chi, parstan_vec)
kpos_donor_ki <- kpos_donor_counts/ndonor_counts

ggplot() +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size=1.3)+
  geom_point(aes(x = data_chi, y = ndonor_counts/(npool_counts * chivec)), size =3, col = 4) + 
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size =2.5) +
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3), labels = ageBMT_names)+
  labs(x = 'Host age', y = NULL, title = paste0("Normalised donor fraction in naive ", toupper(Population), " T cells")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  myTheme + theme(legend.position = c(0.85, 0.15))

ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_counts, col = dataset), size =2.5) + 
  geom_point(aes(x = data_ont, y = ntotal_counts), size =3, col = 4) + 
  geom_point(aes(x = data_chi, y = npool_counts), size =3, col = 4) + 
   scale_color_manual(name= "Age at BMT", values = c(2, 7, 3, 6), labels = ageBMT_names)+
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  #scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + guides(col =F)

ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_kiprop * 100, col = dataset), size =2.5) + 
  geom_point(aes(x = data_ont, y = kpos_pred *100), size =3, col = 4) + 
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3, 6), labels = ageBMT_names)+
  labs(title= paste0("% Ki67high cells in naive ", toupper(Population), " subset"),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_log10(limits = c(0.5, 100))  +
  myTheme + guides(col =F)

kipred_df <- data.frame("donor" = kpos_donor_ki,
                         "host" = kpos_host_ki) %>%
  bind_cols(age.at.S1K = data_chi) %>%
  gather(-c(age.at.S1K), key = subpop, value = ki_prop) 

ggplot() +
  geom_point(data = kipred_df, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size =4) +
  geom_point(data = ki_data, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size =2) + 
  scale_color_manual(name = NULL, values = c(5, 4, 6, 2), labels = c('Donor', "Host")) +
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67high cells in naive ", toupper(Population), " subset")) +
  scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100)) +
  #scale_x_log10(limits = c(5, 450), breaks = c(10, 30, 100, 300, 600)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + theme(legend.position = c(0.85, 0.15))




data_counts %>% filter(ageBMT_bin == "ageBMT_group3") %>%
  summarise(median(age.at.BMT))

chivec1 <- sapply(ts_pred_chi1, Chi_spline)
chivec2 <- sapply(ts_pred_chi2, Chi_spline)
chivec3 <- sapply(ts_pred_chi3, Chi_spline)

## stan solution
ntotal_vec <- N_total_time(ts_pred_ont, parstan_vec)
npool_vec1 <- N_pooled_time(ts_pred_chi1, tb_pred1, parstan_vec)
npool_vec2 <- N_pooled_time(ts_pred_chi2, tb_pred2, parstan_vec)
npool_vec3 <- N_pooled_time(ts_pred_chi3, tb_pred3, parstan_vec)
nhost_vec1 <- N_host_time(ts_pred_chi1, tb_pred1, parstan_vec)
ndonor_vec1 <- N_donor_time(ts_pred_chi1, tb_pred1, parstan_vec)
nhost_vec2 <- N_host_time(ts_pred_chi2, tb_pred2, parstan_vec)
ndonor_vec2 <- N_donor_time(ts_pred_chi2, tb_pred2, parstan_vec)
nhost_vec3 <- N_host_time(ts_pred_chi3, tb_pred3, parstan_vec)
ndonor_vec3 <- N_donor_time(ts_pred_chi3, tb_pred3, parstan_vec)


kpos_vec <- U_total_time(ts_pred_ont, parstan_vec)
kpos_pred <- kpos_vec/ntotal_vec
v1 <- Sys.time()
kpos_pooled_vec <- U_Pooled_time(ts_pred_chi1, tb_pred1, parstan_vec)
v2 <- Sys.time()
v2-v1
kpos_pooled_pred <- kpos_pooled_vec/npool_vec1
kpos_host_vec1 <- U_host_time(ts_pred_chi1, tb_pred1, parstan_vec)
kpos_host_pred1 <- kpos_host_vec1/nhost_vec1
kpos_donor_vec1 <- U_donor_time(ts_pred_chi1, tb_pred1, parstan_vec)
kpos_donor_pred1 <- kpos_donor_vec1/ndonor_vec1
kpos_host_vec2 <- U_host_time(ts_pred_chi2, tb_pred2, parstan_vec)
kpos_host_pred2 <- kpos_host_vec2/nhost_vec2
kpos_donor_vec2 <- U_donor_time(ts_pred_chi2, tb_pred2, parstan_vec)
kpos_donor_pred2 <- kpos_donor_vec2/ndonor_vec2
kpos_host_vec3 <- U_host_time(ts_pred_chi3, tb_pred3, parstan_vec)
kpos_host_pred3 <- kpos_host_vec3/nhost_vec3
kpos_donor_vec3 <- U_donor_time(ts_pred_chi3, tb_pred3, parstan_vec)
kpos_donor_pred3 <- kpos_donor_vec3/ndonor_vec3


p1 <- ggplot() +
  geom_hline(yintercept = 1.0, col = 'darkred', linetype = 2, size=1.3)+
  geom_line(aes(x = ts_pred_chi1, y = ndonor_vec1/(npool_vec1 * chivec1)), size =1.25, col = 2) + 
  geom_line(aes(x = ts_pred_chi2, y = ndonor_vec2/(npool_vec2 * chivec2)), size =1.25, col = 7) + 
  geom_line(aes(x = ts_pred_chi3, y = ndonor_vec3/(npool_vec3 * chivec3)), size =1.25, col = 3) + 
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size =2.5) +
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3), labels = ageBMT_names)+
  labs(x = 'Host age', y = NULL, title = paste0("Normalised donor fraction in naive ", toupper(Population), " T cells")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  myTheme + theme(legend.position = c(0.85, 0.15))

p2 <- ggplot() +
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_counts, col = dataset), size =2.5) + 
  geom_line(aes(x = ts_pred_ont, y = ntotal_vec), size =1.25, col = 6) + 
  geom_line(aes(x = ts_pred_chi1, y = npool_vec1), size =1.25, col = 2) + 
  geom_line(aes(x = ts_pred_chi2, y = npool_vec2), size =1.25, col = 7) + 
  geom_line(aes(x = ts_pred_chi3, y = npool_vec3), size =1.25, col = 3) + 
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3, 6), labels = ageBMT_names)+
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  #scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + guides(col =F)

p3 <- ggplot() +
  geom_line(aes(x = ts_pred_ont, y = kpos_pred * 100), size =1.25, col = 4) + 
  geom_line(aes(x = ts_pred_chi1, y = kpos_pooled_pred * 100), size =1.25, col = 4) + 
  geom_point(data = data_counts, aes(x = age.at.S1K, y = total_kiprop * 100), col =4, size =2) + 
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67high cells in naive ", toupper(Population), " subset")) +
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3, 6), labels = ageBMT_names)+
  scale_y_log10(limits = c(0.1, 100)) +
  scale_x_log10(limits = c(5, 450), breaks = c(10, 30, 100, 300, 600)) +
  myTheme

kidonor_df <- data.frame("ageBMT_group1" = kpos_donor_pred1,
                         "ageBMT_group2" = kpos_donor_pred2,
                         "ageBMT_group3" = kpos_donor_pred3) %>%
  gather(key = ageBMT_bin, value = donor_ki) 

kihost_df <- data.frame("ageBMT_group1" = kpos_host_pred1, 
                        "ageBMT_group2" = kpos_host_pred2, 
                        "ageBMT_group3" = kpos_host_pred3) %>%
  gather(key = ageBMT_bin, value = host_ki)

kipred_df <- data.frame(kidonor_df, kihost_df) %>%
  select(-"ageBMT_bin.1") %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3)) %>%
  gather(-c(age.at.S1K, ageBMT_bin), key = subpop, value = ki_prop) 

p4 <- ggplot() +
  geom_line(data = kipred_df, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size =1.25) +
  geom_point(data = ki_data, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size =2) + 
  scale_color_manual(name = NULL, values = c(2, 4), labels = c('Donor', "Host")) +
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67high cells in naive ", toupper(Population), " subset")) +
  scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100)) +
  #scale_x_log10(limits = c(5, 450), breaks = c(10, 30, 100, 300, 600)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + theme(legend.position = c(0.85, 0.15))


pdf(file = file.path(SaveDir, "PredFig%03d.pdf"),
    width = 12, height = 9, onefile = FALSE, useDingbats = FALSE )

toprow <- cowplot::plot_grid(p2, p1,  labels = c("A", "B"), ncol = 2)
cowplot::plot_grid(toprow, p4,  labels = c("", "C"), nrow = 2)

dev.off()

