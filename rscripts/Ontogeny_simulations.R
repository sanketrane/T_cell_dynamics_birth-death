library(tidyverse)
library(rstan)

#######
Population <-  "cd4"
ModelName <- "asm_deltavar"

WorkinDir <- getwd()
ModelDir <- file.path(WorkinDir, "stan_models/only_chimera", paste0("MAP_", ModelName, "_", Population, ".stan"))
OutputDir <- file.path(WorkinDir, "out_fit/only_chimera", Population, ModelName)
DataDir <- file.path(WorkinDir, "datafiles")



#### importing data for GFP positive proportions
GFPposKineg_df <- readxl::read_excel(file.path(DataDir, "RagGFP_ontogeny_pooled.xlsx"), sheet = 3)%>%
  select('mouseID', 'location', 'age', contains('nai'))
GFPposKipos_df <- readxl::read_excel(file.path(DataDir, "RagGFP_ontogeny_pooled.xlsx"), sheet = 5)%>%
  select('mouseID', 'location', 'age', contains('nai'))


ggplot()+
  geom_point(data = GFPposKineg_df, aes(x=age, y=SP.4nai), col=2)+
  geom_point(data = GFPposKineg_df, aes(x=age, y=LN.4nai), col=4)+
  geom_point(data = GFPposKipos_df, aes(x=age, y=SP.4nai), col=6)+
  geom_point(data = GFPposKipos_df, aes(x=age, y=LN.4nai), col=3)+
  labs(x="Host age", y="GFP %") +ylim(0,100) #+ scale_y_log10(limits=c(0.3, 100))
  

#### Model simulations
rstan::expose_stan_functions(ModelDir)

#Importing fitted parameters
ParamsFile <- read.csv(file.path(OutputDir, paste0("params_", Population, "_", ModelName, ".csv")))
params <- ParamsFile$mean[1:4]
theta <- c(0)
x_i <- c(90)  ## Host age at S1K
x_r <- c(40)  ## Host age at BMT
math_reduce(params, theta, 90, 40)


par_gfp <- c(params, 'r_gfp' = -2.7, 'gfpmax' = 1)

### Vectorise
G_age_psi <- function(age, Time, parms){
  
  mapply(Asm_total_age, age=age, time=Time, MoreArgs = list(parms))
}

Asm_total_age(1, 100, params)
G_age_psi(c(1, 10),  100, par_gfp)

Norm_age_dist <- function(age, Time, parms){
  G_age_psi(age, Time, parms)/
    integrate(G_age_psi, lower = 0, upper = Time, Time=Time, parms=parms)$value
}
Norm_age_dist(1,  100, par_gfp)

normage_df <- data.frame("Host_age_14" = G_age_psi(seq(0, 14, length.out = 100), 14, params),
                         "Host_age_28" = G_age_psi(seq(0, 28, length.out = 100), 28, params),
                         "Host_age_42" = G_age_psi(seq(0, 42, length.out = 100), 42, params),
                         "Host_age_56" = G_age_psi(seq(0, 56, length.out = 100), 56, params),
                         "Host_age_70" = G_age_psi(seq(0, 70, length.out = 100), 70, params))%>%
  gather(key = "Host_age", value = "Counts") %>%
  mutate("cell_age" = c(seq(0, 14, length.out = 100),
                        seq(0, 28, length.out = 100),
                        seq(0, 42, length.out = 100),
                        seq(0, 56, length.out = 100),
                        seq(0, 70, length.out = 100)))

ggplot(normage_df) +
  geom_point(aes(x=cell_age, y=Counts, col=Host_age), size=1.2) +
  labs(x="cell age (days)", title = "cell age distribution at varying host ages")+
  scale_color_discrete(name="Host age", labels = c("2W", "4W", "6W", "8W", "10W"))+
  scale_x_continuous(breaks = c(0, 14, 28, 42, 56, 70))+
  scale_y_log10()
  

#### Deriving the GFP intensity distribution of naive T cells in neonates
G_f_psi <- function(gfp, Time, parms){
  r_gfp = exp(parms[5])
  gfpmax = 1 #exp(parms[6])
  
  ## converting age density into GFP density using a jacobian
  age = (log(gfpmax) - log(gfp))/r_gfp
  value = G_age_psi(age, Time, parms)/(r_gfp * gfp)
  
  return("Gfp_freq" = value)
}

G_f_psi(seq(0, 1, length.out = 100), 10, par_gfp)

Norm_gfp_dist <- function(gfp, Time, parms){
  r_gfp = exp(parms[5])
  gfpmax = 1 #exp(parms[6])
  gfpmin <- exp(-r_gfp * Time)
  
  value = G_f_psi(gfp, Time, parms)/
    integrate(G_f_psi, lower = gfpmin, upper = gfpmax, Time=Time, parms=parms)$value
  
  return(value)
}


normgfp_df <- data.frame("Host_age_14" = G_f_psi(seq(exp(-exp(par_gfp[5])* 14), 1, length.out = 100), 14, par_gfp),
                         "Host_age_28" = G_f_psi(seq(exp(-exp(par_gfp[5])* 28), 1, length.out = 100), 28, par_gfp),
                         "Host_age_42" = G_f_psi(seq(exp(-exp(par_gfp[5])* 42), 1, length.out = 100), 42, par_gfp),
                         "Host_age_56" = G_f_psi(seq(exp(-exp(par_gfp[5])* 56), 1, length.out = 100), 56, par_gfp),
                         "Host_age_70" = G_f_psi(seq(exp(-exp(par_gfp[5])* 70), 1, length.out = 100), 70, par_gfp))%>%
  gather(key = "Host_age", value = "Counts") %>%
  mutate("gfp_intensity" = c(seq(exp(-exp(par_gfp[5])* 14), 1, length.out = 100),
                             seq(exp(-exp(par_gfp[5])* 28), 1, length.out = 100),
                             seq(exp(-exp(par_gfp[5])* 42), 1, length.out = 100),
                             seq(exp(-exp(par_gfp[5])* 56), 1, length.out = 100),
                             seq(exp(-exp(par_gfp[5])* 70), 1, length.out = 100)))

ggplot(normgfp_df) +
  geom_point(aes(x=gfp_intensity, y=Counts, col=Host_age), size=1.2) +
  labs(x="GFP intensity", title = "GFP distribution at varying host ages")+
  scale_color_discrete(name="Host age", labels = c("2W", "4W", "6W", "8W", "10W"))+
  #scale_x_continuous(breaks = c(0, 14, 28, 42, 56, 70))+
  #scale_x_reverse()+
  scale_y_log10()

## Calculating GFP fraction
NofGFPpos <- function(Time, parms){
  gfpmax = 1.0 #exp(parms[6])
  gfp_gate = 0.3
  
  value <- c()
  for (i in 1:length(Time)){
    value[i] <- integrate(G_f_psi, lower = gfp_gate, upper = gfpmax, Time=Time[i], parms=parms)$value
  }
  return(value)
}
NofGFPpos(c(10, 100), par_gfp)

Noftotal <- function(Time, parms){
  gamma_gfp = exp(parms[5])
  gfpmax = 1.0 #exp(parms[6])
  
  value <- c()
  for (i in 1:length(Time)){
    gfpmin <- exp(-gamma_gfp * Time[i])
    value[i] <- integrate(G_f_psi, lower = gfpmin, upper = gfpmax, Time=Time[i], parms=parms)$value
  }
  return(value)
}
Noftotal(c(10, 100),par_gfp)

GFPpos_frac <- function(Time, parms){
  NofGFPpos(Time, parms)/Noftotal(Time, parms)
}

GFPpos_frac(c(1,100), par_gfp)


timeseq <- seq(1, 150, length.out = 100)
GFP_frac <- GFPpos_frac(timeseq, par_gfp)
ggplot()+
  geom_point(aes(x=timeseq, y=GFP_frac), size=1.5)+
  labs(x = "Host age (days)", y="Fraction of GFP positive cells")


## Calculating GFP MFI

### total counts
NofGFPtotal <- function(Time, parms){
  gfpmax = exp(parms[6])
  r_gfp = exp(parms[5])
  
  gfpmin = gfpmax * exp(-r_gfp * Time)
  
  value <- c()
  for (i in 1:length(Time)){
    value[i] <- integrate(G_f_psi, lower = gfpmin, upper = gfpmax, Time=Time[i], parms=parms)$value
  }
  return(value)
}

NofGFPtotal(100, par_gfp)



## Mean GFP intensity
GFP_intensity <- function(gfp, Time, parms){
  gfpmax = exp(parms[6])
  return(gfp * G_f_psi(gfp, Time, parms))
}


GFP_intensity(300, 100, par_gfp)

### Normalised mean GFP intensity
GFP_mfi <- function(Time, parms){
  gfpmax = exp(parms[6])
  r_gfp = exp(parms[5])
  
  # lowest gfp intensity possible at time 'Time'
  gfptime = gfpmax * exp(-r_gfp * Time)
  if(gfptime >= 0.1){
    gfpmin = gfptime
  } else {
    gfpmin = 0.1
  }
  
  answer = integrate(GFP_intensity, lower = gfpmin, upper = gfpmax, Time=Time, parms=parms)$value/
    integrate(G_f_psi, lower = gfpmin, upper = gfpmax, Time=Time, parms=parms)$value
  return(answer)
}


## vectorisation
GFP_mfi_vec <- function(Time, parms){mapply(GFP_mfi, Time, MoreArgs = list(parms))}

GFP_mfi_vec(100, par_gfp)

#### GFP positive cells 
GFP_pos <- function(Time, parms){
  gfpmax = exp(parms[6])
  #gfp_gate = exp(parms[3])
  integrate(GFP_intensity, lower = 1000, upper = gfpmax, Time=Time,  parms=parms)$value/
    integrate(G_f_psi, lower = 1000, upper = gfpmax, Time=Time, parms=parms)$value
}


## vectorisation
GFP_pos_vec <- function(Time, parms){mapply(GFP_pos, Time, MoreArgs = list(parms))}

GFP_pos_vec(100, par_gfp)

