## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
####################################################################################
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[4] = "out.txt"
}

####################################################################################

## model specific details that needs to be change for every run
ModelName <- "asm_rhovar_cd4" #args[1]
Population <- "cd8"# args[2]
hexcode <- args[3]

OutputDir <- file.path("output_csv")
SaveDir <- file.path('out_fit', Population, ModelName) 

## dir to save predictions
if (!file.exists(SaveDir)){
  dir.create(SaveDir)
}

myTheme <- theme(text = element_text(size = 11), axis.text = element_text(size = 11),
                 axis.title =  element_text(size = 11, face = "bold"),
                 plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

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


## dir to save output
if (!file.exists(OutputDir)){
  dir.create(OutputDir)
}
####################################################################################

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file)  %>% 
  arrange(age.at.S1K) %>%
  filter(Nfd <= 1.2)

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34))

data_fit <- chimera_data %>% 
  select(contains('S1K'), contains('counts'), contains('kiprop')) %>%
  mutate(dataset = rep('chimera', nrow(chimera_data))) %>%
  bind_rows(ontogeny_data) %>% arrange(age.at.S1K) 

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file) %>%
  right_join(chimera_data) %>% 
  select(-contains("total"), -contains("Nfd"))%>%
  rename(ageBMTbin=ageBMT_bin)

donor_ki_df <- ki_data %>% filter(subpop == "donor_ki")
host_ki_df <- ki_data %>% filter(subpop == "host_ki")

### reading the models fit # saving output file as 
output_filename=paste(paste0(ModelName,  ".csv"))
model_fit <- read_stan_csv(file.path(OutputDir, "only_chimera", output_filename)) 
## parameters to plot
num_pars <- which(model_fit@model_pars == 'sigma_host_ki')
parameters_to_plot <- model_fit@model_pars[1: num_pars]
nPost = nrow(model_fit)

# ################################################################################################
# calculating PSIS-L00-CV for the fit
## calculating an output from individual runs for the validation of a sucessful run
#loo_loglik1 <- extract_log_lik(model_fit, parameter_name = "log_lik_ont_counts", merge_chains = TRUE)
#loo_loglik2 <- extract_log_lik(model_fit, parameter_name = "log_lik_ont_ki", merge_chains = TRUE)
loo_loglik3 <- extract_log_lik(model_fit, parameter_name = "log_lik_chi_counts", merge_chains = TRUE)
loo_loglik4 <- extract_log_lik(model_fit, parameter_name = "log_lik_Nfd", merge_chains = TRUE)
loo_loglik5 <- extract_log_lik(model_fit, parameter_name = "log_lik_donor_ki", merge_chains = TRUE)
loo_loglik6 <- extract_log_lik(model_fit, parameter_name = "log_lik_host_ki", merge_chains = TRUE)
combined_loglik <- cbind(loo_loglik3, loo_loglik4, loo_loglik5,loo_loglik6)

# optional but recommended
ll_array <- extract_log_lik(model_fit, parameter_name = 'log_lik_chi_counts', merge_chains = FALSE)
r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(combined_loglik,  save_psis = FALSE, cores = 8)

# Widely applicable AIC
AICw_lok <- waic(combined_loglik)

ploocv <- cbind("loo-ic"=loo_loglik$estimates, "WAIC" = AICw_lok$estimates)
ploocv

### parameters table
ptable <- monitor(as.array(model_fit, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
write.csv(out_table, file = file.path(SaveDir, paste0('params_', Population, "_", ModelName, ".csv")))
write.csv(ploocv, file = file.path(SaveDir, paste0('stats_', Population, "_", ModelName, ".csv")))


## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(SaveDir, paste0(Population, "_", ModelName, "QQPlots%03d.pdf")),
    width = 8, height = 5, onefile = FALSE, useDingbats = FALSE)

pairs(model_fit, pars = parameters_to_plot)

### QQ plots
plot(loo_loglik, diagnostic = 'k', label_points = TRUE, main = "PSIS diagnostic plot for k")
plot(loo_loglik, diagnostic = 'n_eff', label_points = TRUE, main = "PSIS diagnostic plot for n_eff")

dev.off()

##### main plots
ageBMT_names <- c(`ageBMT_group1` = "7-9wks",
                  `ageBMT_group2` = "9-11wks", 
                  `ageBMT_group3` = "11-25wks")

ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)

source("rscripts/Chi_tools.R")


pdf(file = file.path(SaveDir, paste0(Population, '_', ModelName, "Fig%03d.pdf")),
    width = 7.2, height = 4.5, onefile = FALSE, useDingbats = FALSE )


#### facet plot for individual age bins
counts_facet <- ggplot() +
  #geom_line(data = Y1pred, aes(x = timeseries, y = median), size =1, col ="#094696") +
  geom_line(data = Y3pred_bin1, aes(x = timeseries, y = median), size =1, col = "#01848a") +
  geom_line(data = Y3pred_bin2, aes(x = timeseries, y = median), size =1, col = 7) +
  geom_line(data = Y3pred_bin3, aes(x = timeseries, y = median), size =1, col = 3) +
  #geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#43b1df", alpha = 0.15) +
  #geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#094696", alpha = 0.3) +
  geom_ribbon(data = Y3pred_bin1, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#01848a", alpha = 0.2) +
  geom_ribbon(data = Y3pred_bin2, aes(x = timeseries, ymin = lb, ymax = ub), fill = 7, alpha = 0.2) +
  geom_ribbon(data = Y3pred_bin3, aes(x = timeseries, ymin = lb, ymax = ub), fill = 3, alpha = 0.2) +
  geom_point(data = ontogeny_data, aes(x = age.at.S1K, y = total_counts), col = "#094696", size=2, alpha = 0.95) +
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = total_counts), col="#01848a",size=2, alpha = 0.95) +
  scale_color_manual(name=NULL, values=c(2, 7, 3)) +
 # scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + guides(col=F)

ki_facet <- ggplot() +
  geom_line(data = Y2pred, aes(x = timeseries, y = median * 100), size =1, col = "#094696") +
  #geom_ribbon(data = kipred, aes(x = timeseries, ymin = lb * 100, ymax = ub * 100), fill = "#43b1df", alpha = 0.15) +
  geom_ribbon(data = Y2pred, aes(x = timeseries, ymin = lb * 100, ymax = ub * 100), fill = "#094696", alpha = 0.3) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_kiprop * 100, col = dataset), alpha = 0.9, size = 2) +
  scale_color_manual(name=NULL, values=c("#01848a","#094696")) +
  labs(title= paste0('% of Ki67+ cell in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(0.3, 100), trans = 'log10') +
  myTheme + theme(legend.position = c(0.85,0.85))

Nfd_plot <- ggplot() +
  geom_hline(yintercept = 1.0, col = "darkred", linetype = 2, size=1.3) +
  geom_line(data = Y4pred_bin1, aes(x = timeseries, y = median), size =1, col = 2) + 
  geom_line(data = Y4pred_bin2, aes(x = timeseries, y = median), size =1, col = 7) + 
  geom_line(data = Y4pred_bin3, aes(x = timeseries, y = median), size =1, col = 3) +
  geom_ribbon(data = Y4pred_bin1, aes(x = timeseries, ymin = lb, ymax = ub), fill = 2, alpha = 0.3) +
  geom_ribbon(data = Y4pred_bin2, aes(x = timeseries, ymin = lb, ymax = ub), fill = 7, alpha = 0.3) +
  geom_ribbon(data = Y4pred_bin3, aes(x = timeseries, ymin = lb, ymax = ub), fill = 3, alpha = 0.3) +
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin), size =2.5) + 
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 3), labels = ageBMT_names)+
  labs(x = 'Host age', y = NULL, title = paste0("Normalised donor fraction in naive ", toupper(Population), " T cells")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  #scale_x_log10(limits = c(50, 450), breaks = c(10, 30, 100, 300, 450)) +
  myTheme + theme(legend.position = c(0.85, 0.15))

ki_facets <- ggplot() +
  geom_line(data = kipred_median, aes(x = age.at.S1K, y = ki_prop*100, col = Popln), size=1.2) +
  geom_ribbon(data = ki_pred_bb, aes(x = age.at.S1K, ymin = lb*100, ymax = ub*100, fill = Popln), alpha = 0.3)+
  geom_point(data = ki_data, aes(x = age.at.S1K, y = ki_prop * 100, col = subpop), size =2) + 
  scale_color_manual(name = NULL, values = c(6, 4), labels = c('Donor', "Host")) +
  scale_fill_manual(name = NULL, values = c(6, 4), labels = c('Donor', "Host")) +
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67high cells in naive ", toupper(Population), " subset")) +
  scale_y_log10(limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100)) +
  #scale_x_log10(limits = c(5, 450), breaks = c(10, 30, 100, 300, 600)) +
  facet_grid(. ~ ageBMTbin, labeller = as_labeller(ageBMT_names)) +
  myTheme + theme(legend.position = c(0.95, 0.9)) 

cowplot::plot_grid(counts_facet, ki_facet, labels = c("A", "B"),  nrow =  1)

toprow <- cowplot::plot_grid(counts_facet,Nfd_plot, labels = c("A", "B"), ncol = 2)
cowplot::plot_grid(toprow, ki_facets,  labels = c("", "C"), nrow = 2)


dev.off()




