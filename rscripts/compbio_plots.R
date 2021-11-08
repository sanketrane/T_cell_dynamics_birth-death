## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
####################################################################################
Population = 'cd8'
## import
DataFile <- file.path("Datafiles", paste0(Population, "_ln.csv"))
FitDir <- file.path("output_csv")
OutputDir <- file.path('out_fit', Population) 

## dir to save output
if (!file.exists(FitDir)){
  dir.create(FitDir)
}
if (!file.exists(OutputDir)){
  dir.create(OutputDir)
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

log10minorbreaks = as.numeric(1:10 %o% 10^(-4:8))

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) 

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

## time course for predictions
data_time <- data_fit$age.at.S1K                        # timepoints in observations 
ts_pred <- seq(min(data_time), max(data_time), length.out = 300)



####################################################################################

counts_cd4 <- ggplot() +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_counts, col = dataset), size = 2) +
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Mouse age (days)") + 
  #scale_x_continuous(limits = c(5, 300), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme 

ki_cd4 <- ggplot() +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_kiprop * 100), col = "darkblue", size = 3) +
  labs(title= paste0('% Ki67+ of cells in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Mouse age (days)") + 
  #scale_x_continuous(limits = c(5, 300), breaks = c(3,10,30,100,300), trans = "log10") +
  scale_y_continuous(limits = c(0.1, 100), trans = 'log10') +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

counts_cd8 <- ggplot() +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_counts), col = "darkblue", size = 2) +
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Mouse age (days)") + 
  #scale_x_continuous(limits = c(5, 300), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 3e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

ki_cd8 <- ggplot() +
   geom_point(data = data_fit, aes(x = age.at.S1K, y = total_kiprop * 100), col = "darkblue", size = 3) +
  labs(title= paste0('% Ki67+ of cells in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Mouse age (days)") + 
  #scale_x_continuous(limits = c(5, 300), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1, 100), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")


cowplot::plot_grid(counts_cd4, counts_cd8, nrow =  1)
cowplot::plot_grid(ki_cd4, ki_cd8,  nrow =  1)


#################################################################
##### neutral model plots
### reading the models fit # saving output file as 
model1_filename = paste0('lipS52', "_", Population, ".csv")
model1_fit <- read_stan_csv(file.path(FitDir, model1_filename))

NY1pred <- as.data.frame(model1_fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

NY2pred <- as.data.frame(model1_fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)


##### lip model plots
model2_filename=paste0('rtemS52', "_", Population, ".csv")
model2_fit <- read_stan_csv(file.path(FitDir, model2_filename))

LY1pred <- as.data.frame(model2_fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

LY2pred <- as.data.frame(model2_fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

##### DD model plots
model3_filename=paste0('asmS52_rhovar', "_", Population, ".csv")
model3_fit <- read_stan_csv(file.path(FitDir, model3_filename))

DY1pred <- as.data.frame(model3_fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)

DY2pred <- as.data.frame(model3_fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)%>%
  filter(timeseries >=5)


moodel_list <- c(paste0('ASM ', '\u03b4', '-varying'), paste0('ASM ', '\u03c1', '-varying'), "LIP")
moodel_list <- c(paste0('LIP'), paste0('RTE'), "ASM Rho-varying")

counts_median_comb <- data.frame("mod1" = NY1pred, "mod2" = LY1pred, "mod3" = DY1pred) %>%
  select(contains('median')) %>%
  'colnames<-' (moodel_list) %>%
  gather(key = "Modelname", value = "Fit")%>%
  bind_cols("timeseries" = rep(ts_pred, 3))

counts_lb_comb <- data.frame("mod1" = NY1pred, "mod2" = LY1pred, "mod3" = DY1pred) %>%
  select(contains('lb')) %>%
  'colnames<-' (moodel_list) %>%
  gather(key = "Modelname", value = "lb")%>%
  bind_cols("timeseries" = rep(ts_pred, 3))

counts_ub_comb <- data.frame("mod1" = NY1pred, "mod2" = LY1pred, "mod3" = DY1pred) %>%
  select(contains('ub')) %>%
  'colnames<-' (moodel_list) %>%
  gather(key = "Modelname", value = "ub")%>%
  bind_cols("timeseries" = rep(ts_pred, 3))

counts_bounds_comb <- left_join(counts_lb_comb, counts_ub_comb, 
                                by = c('timeseries', "Modelname"))

ki_median_comb <- data.frame("mod1" = NY2pred, "mod2" = LY2pred, "mod3" = DY2pred) %>%
  select(contains('median')) %>%
  'colnames<-' (moodel_list) %>%
  gather(key = "Modelname", value = "Fit")%>%
  bind_cols("timeseries" = rep(ts_pred, 3))

ki_lb_comb <- data.frame("mod1" = NY2pred, "mod2" = LY2pred, "mod3" = DY2pred) %>%
  select(contains('lb')) %>%
  'colnames<-' (moodel_list) %>%
  gather(key = "Modelname", value = "lb")%>%
  bind_cols("timeseries" = rep(ts_pred, 3))

ki_ub_comb <- data.frame("mod1" = NY2pred, "mod2" = LY2pred, "mod3" = DY2pred) %>%
  select(contains('ub')) %>%
  'colnames<-' (moodel_list) %>%
  gather(key = "Modelname", value = "ub")%>%
  bind_cols("timeseries" = rep(ts_pred, 3))

ki_bounds_comb <- left_join(ki_lb_comb, ki_ub_comb, 
                                by = c('timeseries', "Modelname"))
head(ki_bounds_comb)

####################################################################################
## fits plots


#### facet plot for individual age bins
counts_facet <- ggplot() +
  geom_line(data = counts_median_comb, aes(x = timeseries, y = Fit, col=Modelname), alpha = 0.8, size =1.2) +
  geom_ribbon(data = counts_bounds_comb, aes(x = timeseries, ymin = lb, ymax = ub, fill = Modelname)
              , alpha = 0.3) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_counts), col = "#A45ee5", alpha = 0.9, size = 2) +
  scale_color_manual(name=NULL, values=c(2,7,3)) +
  scale_fill_manual(name=NULL, values=c(2,7,3))+
  labs(title= paste0('Total numbers of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(1e5, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + theme(legend.position = c(0.85, 0.15))


ki_facet <- ggplot() +
  geom_line(data = ki_median_comb, aes(x = timeseries, y = Fit* 100, col=Modelname), size =1.2) +
  geom_ribbon(data = ki_bounds_comb, aes(x = timeseries, ymin = lb* 100, ymax = ub* 100, fill = Modelname), alpha = 0.2) +
  geom_point(data = data_fit, aes(x = age.at.S1K, y = total_kiprop * 100), alpha = 0.9, col = "#A45ee5", size = 2) +
  scale_color_manual(name=NULL, values=c(2, 7, 3), breaks= moodel_list) +
  scale_fill_manual(name=NULL, values=c(2, 7, 3), breaks= moodel_list)+
  labs(title= paste0('% of Ki67+ cell in naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(5, 450), breaks = c(3,10,30,100,300), trans = "log10")+
  scale_y_continuous(limits = c(0.1, 100), trans = 'log10') +
  myTheme + guides(col=F, fill=F)

cowplot::plot_grid(counts_facet, ki_facet, nrow =  1)


####################################################################################
## params compare
## param dist from the fit
matrix_of_draws1 <- as.data.frame(model1_fit)
matrix_of_draws2 <- as.data.frame(model2_fit)
matrix_of_draws3 <- as.data.frame(model3_fit)

## model1 -- asm_deltavar
alpha1 <- quantile(matrix_of_draws2$N0/4.3e5, probs = c(0.5, 0.025, 0.975))
delta_inv1 <- quantile(matrix_of_draws1$delta, probs = c(0.5, 0.025, 0.975))
rho_inv1 <- quantile(matrix_of_draws1$rho, probs = c(0.5, 0.025, 0.975))
r_del <- quantile(matrix_of_draws1$r_del, probs = c(0.5, 0.025, 0.975))

## model1 -- asm_rhovar
alpha2 <- quantile(matrix_of_draws2$N0/4.3e5, probs = c(0.5, 0.025, 0.975))
delta_inv2 <- quantile(matrix_of_draws2$delta, probs = c(0.5, 0.025, 0.975))
rho_inv2 <- quantile(matrix_of_draws2$rho, probs = c(0.5, 0.025, 0.975))
r_rho <- quantile(matrix_of_draws2$r_rho, probs = c(0.5, 0.025, 0.975))

## model3 -- LIP
alpha3 <- quantile(matrix_of_draws3$psi, probs = c(0.5, 0.025, 0.975))
delta_inv3 <- quantile(matrix_of_draws3$delta, probs = c(0.5, 0.025, 0.975))
rho_inv3 <- quantile(matrix_of_draws3$rho, probs = c(0.5, 0.025, 0.975))

modelnames <- c('asm_del', 'asm_rho', 'ddl')

df_alpha <- data.frame(t(data.frame(alpha1, alpha2, alpha3))) %>%
  mutate(modelName = modelnames,
         par = rep(paste0('Rate of influx','(\u03b1)') , 3))

df_delta <- data.frame(t(data.frame(delta_inv1 , delta_inv2 , delta_inv3))) %>%
  mutate(modelName = modelnames,
         par = rep(paste0('Loss rate','(\u03b4)'), 3))

df_rho <- data.frame(t(data.frame(rho_inv1 , rho_inv2, rho_inv3))) %>%
  mutate(modelName = modelnames,
         par = rep(paste0('Rate of division','(\u03c1)'), 3))



rho_inv3.1 <- quantile(matrix_of_draws3$rho/(1 + (matrix_of_draws3$`y1_mean[3]`/matrix_of_draws3$count_bar)^3), 
                       probs = c(0.5, 0.025, 0.975))
rho_inv3.2 <- quantile(matrix_of_draws3$rho/(1 + (matrix_of_draws3$`y1_mean[22]`/matrix_of_draws3$count_bar)^3), 
                       probs = c(0.5, 0.025, 0.975))
rho_inv3.3 <- quantile(matrix_of_draws3$rho/(1 + (matrix_of_draws3$`y1_mean[34]`/matrix_of_draws3$count_bar)^3), 
                       probs = c(0.5, 0.025, 0.975))


allpars_df <- rbind(df_alpha, df_delta, df_rho)
names(allpars_df) <- c('Estimates', 'par_lb', 'par_ub', 'Model', 'param')

real_modelname <- c(paste0('ASM ', '\u03b4', '-varying'),
                    paste0('ASM ', '\u03c1', '-varying'),
                    "DDL")
                
ggplot(allpars_df)+
  geom_point(aes(y=Estimates, x=Model, col=Model), size = 3) +
  geom_linerange(aes(y=Estimates, ymin=par_lb, ymax=par_ub, x=Model, col=Model),
                 size=0.8,  position=position_dodge(0.05)) +
  scale_color_manual(values = c(4, 3, 2, 2, 2), label = real_modelname)+
  facet_wrap(~ param) + scale_y_log10() + 
  myTheme + theme(axis.text.x=element_blank(),
                  axis.title.x=element_blank())
  
  #labs(y='Parameters', x= "Estimates") +
  #scale_y_continuous(limits = c(1e-4, 1.2), trans="log10", breaks=c(3e-3, 3e-2, 3e-1),
                 #    minor_breaks = log10minorbreaks, labels =fancy_scientific) +



#pars_to_plot <- c('\u03b1', '\u03b4', '\u03c1', 'r-\u03c1')

#plot_table <- data.frame(t(data.frame(alpha1, delta1, rho1))) %>%
# mutate(pars_plot = pars_to_plot)
plot_table <- data.frame(t(data.frame(alpha1, delta1, rho1, r_del1))) %>%
  mutate(pars_plot = pars_to_plot)

names(plot_table) <- c('pars_median', 'pars_lb', 'pars_ub', 'pars_plot')

ggplot(plot_table)+
  geom_point(aes(x=pars_median, y=pars_plot), size = 3)+
  geom_linerange(aes(x=pars_median, xmin=pars_lb, xmax=pars_ub, y=pars_plot),
                 size=0.8,  position=position_dodge(0.05)) +
  labs(y='Parameters', x= "Estimates") +
  scale_x_continuous(limits = c(1e-3, 1.0), trans="log10", breaks=c(3e-3, 3e-2, 3e-1),
                     minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme


ggplot(out_table2)+
  geom_point(aes(x=pars_median, y=parameters_to_plot), size = 3)+
  geom_linerange(aes(x=pars_median, xmin=pars_lb, xmax=pars_ub, y=parameters_to_plot),
                 size=0.8,  position=position_dodge(0.05)) +
  scale_x_log10(limits=c(0.001, 3))

write.csv(out_table, file = file.path(TablDir, paste0('params_', Population, "_", modelName, ".csv")))

precs_m1 <- rethinking::precis(model_fit, pars = parameters_to_plot[2:5])
plot(out_table)




