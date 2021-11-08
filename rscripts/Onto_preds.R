library(tidyverse)
library(rstan)
library(parallel)

#######
Population <-  'cd4'
ModelName <- "asm_deltavar_cd4_K50"

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) %>% 
  arrange(age.at.S1K) %>%
  filter(Nfd <= 1.2)

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34)) %>% 
  arrange(age.at.S1K) 



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



#################################

# Extracting parameters from the model fit
### reading the models fit # saving output file as 
output_filename=paste(paste0(ModelName, "_", Population, ".csv"))
model_fit <- read_stan_csv(file.path("output_csv", "only_chimera", paste0(ModelName, ".csv"))) 
## parameters to plot
num_pars <- which(model_fit@model_pars == 'r_del')
parameters_to_plot <- model_fit@model_pars[1: num_pars]
nPost = nrow(model_fit)

# dataframe with each row as samoled during the fit
matrix_of_draws <- as.data.frame(model_fit, pars = parameters_to_plot)

# stan model
rstan::expose_stan_functions(file.path("stan_models", "only_onto", paste0("MAP_", ModelName, ".stan")))

# time points for predictions
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 100)

# predictions
tim_st <- Sys.time()
total_counts <- data.frame()
total_ki <- data.frame()
num_cores <- detectCores() - 6
for (i in 1:nPost){
  params <- unlist(matrix_of_draws[i, ])
  total_counts[i, 1:100] <- N_total_time(ts_pred_ont, params)
  total_ki[i, 1:100] <- mcmapply(U_total_time,  ts_pred_ont, MoreArgs =list(params), mc.cores = num_cores)/total_counts[i, 1:100] 
}
tim_end <- Sys.time()
tot_time <- tim_end - tim_st
print(tot_time)

Y1pred <- total_counts %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_ont)%>%
  filter(timeseries >=5)

Y2pred <- total_ki %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_ont)%>%
  filter(timeseries >=5)

print("DONE!")

Pred_dir <- file.path('out_fit', 'Onto_preds', Population, ModelName) 

## dir to save predictions
if (!file.exists(Pred_dir)){
  dir.create(Pred_dir)
}

pdf(file = file.path(Pred_dir, paste0(Population, '_', ModelName, "weFig%03d.pdf")),
    width = 12, height = 4.5, onefile = FALSE, useDingbats = FALSE )


counts_facet <- ggplot() +
  geom_line(data = Y1pred, aes(x = timeseries, y = median), size =1, col ="#094696") +
  geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#094696", alpha = 0.3) +
  geom_point(data = ontogeny_data, aes(x = age.at.S1K, y = total_counts), col = "#094696", size=2, alpha = 0.95) +
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = total_counts), col = 2, size=2, alpha = 0.95) +
  scale_color_manual(name=NULL, values=c(2, 7, 3)) +
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 5e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + guides(col=F)

ki_facet <- ggplot() +
  geom_line(data = Y2pred, aes(x = timeseries, y = median * 100), size =1, col = "#094696") +
  geom_ribbon(data = Y2pred, aes(x = timeseries, ymin = lb * 100, ymax = ub * 100), fill = "#094696", alpha = 0.3) +
  geom_point(data = ontogeny_data, aes(x = age.at.S1K, y = total_kiprop * 100), col = "#094696", alpha = 0.9, size = 2) +
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = total_kiprop * 100), col = 2, alpha = 0.9, size = 2) +
  scale_color_manual(name=NULL, values=c("#01848a","#094696")) +
  labs(title= paste0('% of Ki67+ cell in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(0.3, 100), trans = 'log10') +
  myTheme + theme(legend.position = c(0.85,0.85))

cowplot::plot_grid(counts_facet, ki_facet, align = c("h"))
dev.off()
