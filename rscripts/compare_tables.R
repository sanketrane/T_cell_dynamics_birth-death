rm(list = ls())  
gc()    

library(tidyverse)
library(gridExtra)
library(formattable)


####################################################################################
Population = 'cd8'
## directories for saving outputs
OutputDir <- file.path('out_fit', "only_chimera", Population) 
FigDir <- file.path(OutputDir, 'figures')
TablDir <- file.path(OutputDir, 'tables')
LooDir <- file.path('loo_fit', "only_chimera") 

## dir to save output
if (!file.exists(OutputDir)){
  dir.create(OutputDir)
}
if (!file.exists(FigDir)){
  dir.create(FigDir)
}
if (!file.exists(TablDir)){
  dir.create(TablDir)
}

customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customOrange = "#FF8C00"
customRed = "#FF6347"

## function for exporting a data table of delta loo ic values and akaike weights
#takes 2 separate lists of the name of the models and loo-ic values
model_compare <- function(looiclist){
  # delta loo-ic
  deltaloo_list <- looiclist - min(looiclist)
  ## akaike wt
  akaikewt_numerator <- function(x) exp(-0.5 * x)
  akaikewt_list <- sapply(deltaloo_list, akaikewt_numerator) * 100/sum(exp(- 0.5 * deltaloo_list))
  
  export_table <- data.frame('deltaloo' = round(deltaloo_list, 2),
                             'Akaike_wt' = round(akaikewt_list, 2))
  colnames(export_table)[1:2] <-  c(paste0('\u0394', 'LooIC'),  paste0('Akaike weight ', '\u0025'))
  
  return(export_table)
}

### reading the loo objects for each model
neutral_filename = paste0('loosave_', 'neutral', '_', Population,    '.rds')
lip_filename = paste0('loosave_',  'lip', '_', Population,   '.rds')
ddm_filename = paste0('loosave_', 'ddm', '_', Population,  '.rds')
rtem_filename = paste0('loosave_', 'rtem', '_',Population,    '.rds')
rtemld_filename = paste0('loosave_',  'rtemld', '_', Population,  '.rds')
asm_deltavar_filename = paste0('loosave_', 'asm_deltavar', '_', Population,    '.rds')
asm_rhovar_filename = paste0('loosave_', 'asm_rhovar', '_', Population, '.rds')


neutral_loo <- readRDS(file.path(LooDir, neutral_filename))
lip_loo <- readRDS(file.path(LooDir, lip_filename))
ddm_loo <- readRDS(file.path(LooDir, ddm_filename))
rtem_loo <- readRDS(file.path(LooDir, rtem_filename))
rtemld_loo <- readRDS(file.path(LooDir, rtemld_filename))
asm_deltavar_loo <- readRDS(file.path(LooDir, asm_deltavar_filename))
asm_rhovar_loo <- readRDS(file.path(LooDir, asm_rhovar_filename))

model_list <- list('Neutral' = neutral_loo, 
                   'LIP' = lip_loo, 
                   "DD" = ddm_loo, 
                   "RTE" = rtem_loo,
                   "RTE-MLD" = rtemld_loo,
                   "ASM_deltavar" = asm_deltavar_loo,
                   "ASM_rhovar" = asm_rhovar_loo)
comp <- loo_compare(model_list)
print(comp, simplify = F)
elpd_list <- comp[, 7] - min(comp[, 7])
akaikewt_numerator <- function(x) exp(-0.5 * x)
akaikewt_list <- sapply(elpd_list, akaikewt_numerator) * 100/ sum(exp(- 0.5 *elpd_list))


export_table <- formattable(model_compare(elpd_list),
                             align = c('l', 'l', 'c', 'c'),
                             list(`Akaike weight %` =
                                    formatter("span", style = ~ style(font.size='16px')),
                                  `ΔLooIC` =
                                    formatter("span", style = x ~ ifelse(x == 0.00,
                                                                         style(color = customGreen, font.size='16px'),
                                                                         ifelse(x <= 5.0, style(color = customOrange, font.size='16px'), 
                                                                                ifelse(x <100, style(font.size='16px'), NA)))))
)

export_table
loo_model_weights(model_list, method = 'pseudobma') * 100


######################################################################################################
######################################################################################################



