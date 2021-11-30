#clearing R environment 
rm(list = ls()); gc()

# imorting libraries
library(tidyverse)
library(readxl)
library(parallel)


## model and data specific details
Population <-  "cd8"
ModelName <- "asm_deltavar"

## setting directories to work with 
#setwd("~/gfp_preds")
WorkinDir <- getwd()
DataDir <- file.path(WorkinDir, "datafiles")
WorkinDir <- getwd()
DataDir <- file.path(WorkinDir, "datafiles")
ModelDir <- file.path(WorkinDir, "stan_models", paste0("MAP_", ModelName, "_", Population, ".stan"))
OutputDir <- file.path(WorkinDir, Population, ModelName)

### Theme for plots
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), 
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.text = element_text(size=12), legend.title = element_text(12))

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



#### importing data
# GFP positive fraction within Ki67 negative population
GFPposKineg_df <- readxl::read_excel(file.path(DataDir, "RagGFP_ontogeny_pooled2.xlsx"), sheet = 3)%>%
  select('mouseID', 'location', 'age', 'TH.SP8', "SP.8nai", 'LN.8nai') 

# GFP positive fraction within Ki67 positive population
GFPposKipos_df <- readxl::read_excel(file.path(DataDir, "RagGFP_ontogeny_pooled2.xlsx"), sheet = 6)%>%
  select('mouseID', 'location', 'age', 'TH.SP8',  "SP.8nai", 'LN.8nai') 

# Total GFP fraction
GFP_pos_total <- full_join(GFPposKineg_df, GFPposKipos_df, 
                           by=c("mouseID", "location", "age")) %>%
  mutate(gfp_pos_Thy = TH.SP8.x + TH.SP8.y,
         gfp_pos_LN = LN.8nai.x + LN.8nai.y,
         gfp_pos_Spl = SP.8nai.x + SP.8nai.y)%>%
  select(-contains('.x'), -contains('.y')) %>%
  gather(-c('mouseID', 'location', 'age'), key = 'population', value = 'GFP_percents') %>%
  mutate(GFP_frac = GFP_percents) %>%
  arrange(age)

p0 <- ggplot(GFP_pos_total)+
  geom_point(aes(x=age, y=(GFP_percents), col=population), size=2)+
  labs(x="Host age", y=NULL, title="GFP % total cd8 T cells") + ylim(0, 100)


### Age structured model

#Importing fitted parameters
ParamsFile <- read.csv(file.path(OutputDir, paste0("params_", Population, "_", ModelName, ".csv")))
params_imported <- ParamsFile$mean[1:4]

#### import the age structured model in R from stanfile 
rstan::expose_stan_functions(ModelDir)

### testing the model
# age distribution from the stan model
test1 <- Asm_total_age(age=1, time=100, parms=params_imported)

### Vectorising the stan derived age distrubution to work with in R
G_age_Vec <- function(age, Time, parms){
  mapply(Asm_total_age, age=age, time=Time, MoreArgs = list(parms))
}

test2 <- G_age_Vec(age=1, Time=100, parms=params_imported)
if(test1 != test2){
  print("Error in vectorization of the age distribution")
}

# Normalised age distribution
Norm_age_dist <- function(age, Time, parms){
  G_age_Vec(age, Time, parms)/
    integrate(G_age_Vec, lower = 0, upper = Time, Time=Time, parms=parms)$value
}

normage_df <- data.frame("Host_age_14" = G_age_Vec(seq(0, 14, length.out = 100), 14, params_imported),
                         "Host_age_28" = G_age_Vec(seq(0, 28, length.out = 100), 28, params_imported),
                         "Host_age_42" = G_age_Vec(seq(0, 42, length.out = 100), 42, params_imported),
                         "Host_age_56" = G_age_Vec(seq(0, 56, length.out = 100), 56, params_imported),
                         "Host_age_70" = G_age_Vec(seq(0, 70, length.out = 100), 70, params_imported))%>%
  gather(key = "Host_age", value = "Counts") %>%
  mutate("cell_age" = c(seq(0, 14, length.out = 100),
                        seq(0, 28, length.out = 100),
                        seq(0, 42, length.out = 100),
                        seq(0, 56, length.out = 100),
                        seq(0, 70, length.out = 100)))

p_age <- ggplot(normage_df) +
  geom_point(aes(x=cell_age, y=Counts, col=Host_age), size=1.2) +
  labs(x="cell age (days)", title = "Cell age distribution at varying host ages")+
  scale_color_discrete(name="Host age", labels = c("2W", "4W", "6W", "8W", "10W"))+
  scale_x_continuous(breaks = c(0, 14, 28, 42, 56, 70))

## Calculating GFP fraction
#parameters that control gfp distribution
gfpmax <- 1.0    ### maximum gfp intensity normalized to 1 but can be moified to the observed value.
gamma_gfp <- log(0.067)   ### rate of loss  of gfp expression
gfp_gate <- log(0.3)       ### threshold gfp intensity above which cells considered GFP+
a_bar <- log(10)
par_gfp <- c(a_bar) #c(gamma_gfp, gfp_gate) #,gfpmax)
params <- c(params_imported, par_gfp)

## number of gfp positive cells. gfp intensities >= gfp_gate.
NofGFPpos <- function(Time, parms){
  #gamma_gfp = exp(parms[5])
  #gfp_gate = exp(parms[6])
  #gfpmax = 1.0 #exp(parms[7])
  a_bar = exp(parms[5])
  #a_bar = (1/gamma_gfp) * -log(gfp_gate)
  value <- c()
  for (i in 1:length(Time)){
    value[i] <- ifelse(a_bar <= Time[i],
                       integrate(G_age_Vec, lower = 0, upper = a_bar, Time=Time[i], parms=parms)$value,
                       integrate(G_age_Vec, lower = 0, upper = Time[i], Time=Time[i], parms=parms)$value
    )
  }
  return(value)
}

## Total number of cells within the population
Noftotal <- function(Time, parms){
 
  value <- c()
  for (i in 1:length(Time)){
    gfpmin <- exp(-gamma_gfp * Time[i])
    value[i] <- integrate(G_age_Vec, lower = 0, upper = Time[i], Time=Time[i], parms=parms)$value
  }
  return(value)
}

## fraction of cells that GFP+
GFPpos_frac <- function(Time, parms){
  NofGFPpos(Time, parms)/Noftotal(Time, parms)
}

## data transformation function
logit_transf <- function(p){
  log(p/(1.2-p))
}

## data for fitting
data_fit <- GFP_pos_total %>%
  filter(population == "gfp_pos_LN") %>%
  mutate(subset = "total_pop")%>%
  select('mouseID', 'age', 'population', 'subset', 'GFP_frac') 


GFPfrac_fit <- function(Time, a_bar){
  param <- c(params_imported, a_bar)
  GFPpos_frac(Time, param)
}

nls_gfp <- nls(logit_transf(GFP_frac/100) ~ logit_transf(GFPfrac_fit(age, a_bar)), 
               data = data_fit, start = list('a_bar' = 2.3))
summary(nls_gfp)

par_estm <- c(params_imported, coef(nls_gfp))  #c(params_imported, fit_gfp_frac$par)
## time takes for cell to fall out of the GFP gate in our experiments
a_bar = exp(par_estm[5])
print(paste0("Time taken for a cell with maximum GFP intensity to fall out of the GFP positive gate in our experiments = ", round(a_bar, 2), " days"))
#print(paste0("Estimated GFP half life on naive ", Population, ' is ', round(log(2)/(par_estm[5]), 2), " days"))

## generating predictions from the model fit
Timeseq <- seq(from = 5,
               to = 450,
               length.out=500)

GFP_frac_fit <- GFPpos_frac(Timeseq, par_estm)

p1 <- ggplot()+
  geom_line(aes(x=Timeseq, y=GFP_frac_fit), col=2)+
  geom_point(data=data_fit, aes(x=age, y=GFP_frac), col=3, size=2) +
  ylim(0,1)+ #scale_color_discrete(name=NULL, labels=c("LN", "spleen"))+
  scale_x_log10() + guides(col='none')+
  labs(x = "Host age (days)", y=NULL, title = "ASM fit to the total GFP+fraction in naive CD4 T cells") + myTheme




#### Predicting GFP fractions in Ki67- and Ki67+ compartments
## wrangling the data
GFPposKineg_grouped <- GFPposKineg_df %>%
  select(-'location') %>%
  gather(-c(mouseID, age), key = 'population', value = 'GFP_frac')

GFPposKipos_grouped <- GFPposKipos_df %>%
  select(-'location') %>%
  gather(-c(mouseID, age), key = 'population', value = 'GFP_frac')


GFP_ki_split <- full_join(GFPposKineg_grouped, GFPposKipos_grouped, 
                          by=c("mouseID", "age", "population"),
                          suffix=c('.kineg', '.kipos'))  %>%
  gather(-c(mouseID, age, population), key = 'subset', value = 'GFP_frac') %>%
  filter(population == "LN.8nai") %>%
  rbind(data_fit)


grid_labels <- as_labeller(c('total_pop' = "Total GFP+", 'GFP_frac.kineg' = "GFP+ Ki67-", 'GFP_frac.kipos' = "GFP+ Ki67+"))
p3 <- ggplot(GFP_ki_split)+
  geom_point(aes(x=age, y=(GFP_frac), col=Population), size=2)+
  labs(x="Host age", y=NULL, title="GFP % in CD8 T cells") + 
  #scale_x_log10() +
  #scale_y_log10()+
  facet_wrap(.~ Subset, labeller = grid_labels)

## Proportion of Ki67+ cells 
ki67_pos_frac <- function(Time, parms){
  value = U_total_time(Time, parms)/
    N_total_time(Time, parms)
  return(value)
} 

Counts_pred <- N_total_time(Timeseq, params_imported)
Kipos_frac <- ki67_pos_frac(Timeseq, params_imported)
GFPpos_kipos_preds <- GFPpos_frac(Timeseq, par_estm) * Kipos_frac
GFPpos_kineg_preds <- GFPpos_frac(Timeseq, par_estm) * (1 - Kipos_frac)

GFP_split_pred_df <- data.frame("Host_age" = Timeseq,
                                "total_pop" = GFP_frac_fit,
                                "GFP_frac.kipos" = GFPpos_kipos_preds,
                                "GFP_frac.kineg" = GFPpos_kineg_preds) %>%
  gather(-c(Host_age), key = 'subset', value = 'GFP_frac') %>%
  filter(Host_age <=120) %>%filter(Host_age >=10)


GFP_ki_split$subset <- factor(GFP_ki_split$subset, levels=c("total_pop","GFP_frac.kineg","GFP_frac.kipos")) 
GFP_split_pred_df$subset <- factor(GFP_split_pred_df$subset, levels=c("total_pop","GFP_frac.kineg","GFP_frac.kipos")) 

ggplot()+
  geom_point(data=GFP_ki_split, aes(x=age, y=(GFP_frac)),col="#037971", size=2)+
  labs(x="Host age", y=NULL, title="Predictions of GFP % in CD8 T cells") + 
  geom_line(data=GFP_split_pred_df, aes(x=Host_age, y=GFP_frac*100), col=2)+
  facet_wrap(.~ subset, labeller = grid_labels) +
  scale_x_log10() +myTheme 


#ggsave(filename = file.path(OutputDir, "gfp_fit.pdf"), p1, device = "pdf", width = 4.5, height = 3.5)
ggsave(filename = file.path(OutputDir, "gfp_preds.pdf"), last_plot(), device = "pdf", width = 9, height = 3.0)


## Total Ki67+ proportions 
# GFP positive fraction within Ki67 positive population
GFPnegKipos_df <- readxl::read_excel(file.path(DataDir, "RagGFP_ontogeny_pooled2.xlsx"), sheet = 5)%>%
  select('mouseID', 'location', 'age', 'TH.SP8',  "SP.8nai", 'LN.8nai') 


GFP_numbers_df <- readxl::read_excel(file.path(DataDir, "RagGFP_ontogeny_pooled2.xlsx"), sheet = 1)%>%
  select('mouseID', 'age', 'LN.8nai', 'SP.8nai') %>%
  mutate(age.at.S1K = age,
         total_counts = LN.8nai + SP.8nai) %>%
  select('mouseID', contains('S1K'), contains('counts'))

GFP_data <- full_join(GFPnegKipos_df, GFPposKipos_df, by = c('mouseID', 'location', 'age'))%>%
  mutate(Ki_Thy = (TH.SP8.x + TH.SP8.y)/100,
         Ki_LN= (LN.8nai.x + LN.8nai.y)/100,
         Ki_Spl = (SP.8nai.x + SP.8nai.y)/100) %>%
  select('mouseID', contains('age'), contains('Ki_LN')) %>%
  rename(age.at.S1K = age,
         total_kiprop = Ki_LN)%>%
  left_join(GFP_numbers_df) %>%
  select(-contains('ID')) %>%
  mutate(dataset = rep('RAG_GFP', 19))



## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file)  %>% 
  arrange(age.at.S1K) %>%
  filter(Nfd <= 1.2) %>%
  select(contains('S1K'), contains('counts'), contains('kiprop'))

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  select(contains('S1K'), contains('counts'), contains('kiprop'))%>%
  mutate(dataset = rep('ontogeny', 34))

# reading and reshaping data -------------------------------------------------
Onto_counts <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 2)
Onto_ki67 <- read_excel("datafiles/ontogeny_data.xlsx", sheet = 3)

counts_ln  <- Onto_counts %>%
  filter(organ == "LN") %>%
  select(mouse.id, age.at.S1K.days, CD8.nai)%>%
  rename(age.at.S1K = age.at.S1K.days, 
         counts = CD8.nai)

ki67_ln  <- Onto_ki67 %>%
  filter(organ == "LN") %>%
  select(mouse.id, age.at.S1K.days, CD8.nai)%>%
  rename(age.at.S1K = age.at.S1K.days, 
         ki67 = CD8.nai)

counts_spl  <- Onto_counts %>%
  filter(organ == "SP") %>%
  select(mouse.id, age.at.S1K.days, CD8.nai) %>%
  rename(age.at.S1K = age.at.S1K.days, 
         counts = CD8.nai)

Onto_data <- full_join(counts_spl, counts_ln, by=c('mouse.id', 'age.at.S1K'), suffix=c('.splthrow', '.lnthrow')) %>%
  na.omit() %>%
  mutate(total_counts = counts.lnthrow + counts.splthrow) %>%
  left_join(ki67_ln, by=c('mouse.id', 'age.at.S1K'))%>%
  mutate(dataset = rep('ontogeny', 30),
         total_kiprop = ki67/100) %>%
  select(-contains('throw'), -mouse.id, -ki67)


data_pred <- chimera_data %>% 
  mutate(dataset = rep('chimera', nrow(chimera_data))) %>%
  bind_rows(Onto_data) %>% 
  bind_rows(GFP_data) %>% 
  arrange(age.at.S1K) 

ggplot()+
  labs(x="Host age", y=NULL, title="Ki67 % in CD8 T cells") + 
  scale_x_log10()+ scale_y_log10()+
  scale_color_manual(name=NULL, values=c(2,"#0571c4","#037971"))+
  geom_line(aes(x=Timeseq, y=Kipos_frac*100), col="2", size=0.7) + 
  geom_point(data=data_pred, aes(x=age.at.S1K, y=(total_kiprop * 100), col=dataset), size=2)+
  myTheme + guides(col='none')

ggsave(filename = file.path(OutputDir, "Ki67_pred.pdf"), last_plot(), device = "pdf", width = 4.8, height = 3.45)


ggplot()+
  labs(x="Host age", y=NULL, title="Total counts of CD8 T cells") + 
  scale_x_log10()+ 
  scale_color_manual(name=NULL, values=c(2,"#0571c4","#037971"))+
  scale_y_continuous(limits = c(4e5, 5e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  geom_line(aes(x=Timeseq, y=Counts_pred), col=2, size=0.7) +
  geom_point(data=data_pred, aes(x=age.at.S1K, y=(total_counts), col=dataset), size=2)+
  myTheme + guides(col='none')

ggsave(filename = file.path(OutputDir, "counts_pred.pdf"), last_plot(), device = "pdf", width = 4.8, height = 3.45)


##Host age effects
hostageDir <- file.path(WorkinDir, "stan_models", paste0("asm_da_dh_cd8.stan"))
rstan::expose_stan_functions(hostageDir)

params_HA <- c(params_imported, 15)

lambda_h <- lam
ggplot()+
  geom_line(aes(x=Timeseq-4, y=delta_time(Timeseq, params_HA)))+
  geom_line(aes(x=Timeseq-4, y=delta_time(Timeseq, c(params_imported, 15, 4))), col=2)+
  labs(x="Host age", y=NULL, title="Rate of turnover of CD8 RTE (delta0)") + 
  scale_x_log10(breaks=c(1, 3, 10, 30, 100, 300)) + myTheme

ggsave(filename = file.path(OutputDir, "loss_rate_RTE.pdf"), last_plot(), device = "pdf", width = 3.0, height = 3.0)


Counts_HA_frac <- N_total_time(Timeseq, params_HA)
Kipos_HA_frac <- unlist(mclapply(Timeseq, ki67_pos_frac, parms=params_HA, mc.cores = 4))


GFPpos_kipos_HA_preds <- GFPpos_frac(Timeseq, par_estm) * Kipos_HA_frac
GFPpos_kineg_HA_preds <- GFPpos_frac(Timeseq, par_estm) * (1 - Kipos_HA_frac)

GFP_HA_pred_df <- data.frame("Host_age" = Timeseq,
                             "GFP_frac.kipos" = GFPpos_kipos_HA_preds,
                             "GFP_frac.kineg" = GFPpos_kineg_HA_preds) %>%
  gather(-c(Host_age), key = 'subset', value = 'GFP_frac')%>%
  filter(Host_age <=120) %>%filter(Host_age >=10)

GFP_split_pred_df$subset <- factor(GFP_split_pred_df$subset, levels=c("GFP_frac.kineg","GFP_frac.kipos")) 


##Host age effects
hostageDir <- file.path(WorkinDir, "stan_models", paste0("asm_da_rh_cd8.stan"))
rstan::expose_stan_functions(hostageDir)

params_HA2 <- c(params_imported, 15)


Counts_HA_frac2 <- N_total_time(Timeseq, params_HA)
Kipos_HA_frac2 <- unlist(mclapply(Timeseq, ki67_pos_frac, parms=params_HA, mc.cores = 4))

ggplot()+
  labs(x="Host age", y=NULL, title="Ki67 % in CD8 T cells") + 
  scale_x_log10()+ scale_y_log10()+
  scale_color_manual(name=NULL, values=c(2,"#0571c4","#037971")) +
  geom_line(aes(x=Timeseq, y=Kipos_frac*100), col=2) +
  geom_line(aes(x=Timeseq, y=Kipos_HA_frac*100), col="#ddaadd", size=2) +
  geom_line(aes(x=Timeseq, y=Kipos_HA_frac2*100), col=4) +
  geom_point(data=data_pred, aes(x=age.at.S1K, y=(total_kiprop * 100), col=dataset), size=1.52)+
  myTheme + guides(col='none')

#ggsave(filename = file.path(OutputDir, "Ki67pred_Hostage.pdf"), last_plot(), device = "pdf", width = 4.8, height = 3.45)


ggplot()+
  labs(x="Host age", y=NULL, title="Total counts of CD8 T cells") + 
  scale_x_log10()+ 
  scale_color_manual(name=NULL, values=c(2,"#0571c4","#037971"))+
  scale_y_continuous(limits = c(3e5, 5e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  geom_line(aes(x=Timeseq, y=Counts_pred), col=2) +
  geom_line(aes(x=Timeseq, y=Counts_HA_frac), col="#ddaadd", size=2) +
  geom_line(aes(x=Timeseq, y=Counts_HA_frac2), col=4) +
  geom_point(data=data_pred, aes(x=age.at.S1K, y=total_counts, col=dataset), size=1.52)+
  myTheme + guides(col='none')


GFPpos_kipos_HA_preds2 <- GFPpos_frac(Timeseq, par_estm) * Kipos_HA_frac2
GFPpos_kineg_HA_preds2 <- GFPpos_frac(Timeseq, par_estm) * (1 - Kipos_HA_frac2)

GFP_HA_pred_df2 <- data.frame("Host_age" = Timeseq,
                             "GFP_frac.kipos" = GFPpos_kipos_HA_preds2,
                             "GFP_frac.kineg" = GFPpos_kineg_HA_preds2) %>%
  gather(-c(Host_age), key = 'subset', value = 'GFP_frac')%>%
  filter(Host_age <=120) %>%filter(Host_age >=10)

GFP_split_pred_df$subset <- factor(GFP_split_pred_df$subset, levels=c("GFP_frac.kineg","GFP_frac.kipos")) 


ggplot()+
  geom_point(data=filter(GFP_ki_split, subset!='total_pop'), aes(x=age, y=(GFP_frac)),col="#037971", size=2)+
  labs(x="Host age", y=NULL, title="Predictions of GFP % in CD8 T cells") + 
  geom_line(data=filter(GFP_split_pred_df, subset!='total_pop'), aes(x=Host_age, y=(GFP_frac*100)), col=2)+
  geom_line(data=GFP_HA_pred_df, aes(x=Host_age, y=GFP_frac*100), col="#ddaadd", size=2)+
  geom_line(data=GFP_HA_pred_df2, aes(x=Host_age, y=GFP_frac*100), col=4)+
  facet_wrap(.~ subset, labeller = grid_labels) +
  scale_x_log10()  +myTheme

ggsave(filename = file.path(OutputDir, "gfppreds_Hostage.pdf"), last_plot(), device = "pdf", width = 6, height = 3.0)

