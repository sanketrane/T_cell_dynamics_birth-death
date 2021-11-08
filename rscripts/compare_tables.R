rm(list = ls())  
gc()    

library(tidyverse)
library(gridExtra)
library(formattable)

## directories for saving outputs
OutputDir <- file.path('compbio_talk')
FigDir <- file.path(OutputDir, 'figures')
TablDir <- file.path(OutputDir, 'tables')

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
neutral_filename = paste0('loo_', Population, '_', 'neutral',  '.rds')
lip_filename = paste0('loo_', Population, '_', 'lip',  '.rds')
ddm_filename = paste0('loo_', Population, '_', 'ddm',  '.rds')
rtem_filename = paste0('loo_', Population, '_', 'rtem',  '.rds')
rtemld_filename = paste0('loo_', Population, '_', 'rtemld',  '.rds')
asm_deltavar_filename = paste0('loo_', Population, '_', 'asm_deltavar',  '.rds')
asm_rhovar_filename = paste0('loo_', Population, '_', 'asm_rhovar',  '.rds')


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



elpd_list <- comp[, 7] - min(comp[, 7])
akaikewt_numerator <- function(x) exp(-0.5 * x)
akaikewt_list <- sapply(elpd_list, akaikewt_numerator) * 100/ sum(exp(- 0.5 *elpd_list))

loo_model_weights(model_list, method = 'pseudobma') * 100


######################################################################################################
######################################################################################################

model_compare <- function(modellist, looiclist){
  # delta loo-ic
  deltaloo_list <- looiclist - min(looiclist)
  ## akaike wt
  akaikewt_numerator <- function(x) exp(-0.5 * x)
  akaikewt_list <- sapply(deltaloo_list, akaikewt_numerator) * 100/sum(exp(- 0.5 * deltaloo_list))
  
  export_table <- data.frame('Model' = modellist,
                             'deltaloo' = round(deltaloo_list, 2),
                             'Akaike_wt' = round(akaikewt_list, 2))
  colnames(export_table)[2:3] <-  c(paste0('\u0394', 'LooIC'),  paste0('Akaike weight ', '\u0025'))
  
  return(export_table)
}


models_list2 <- c('Neutral', 'LIP', 'Density-dependent loss', 'RTE', 'RTE-MLD', "ASM Delta varying", "ASM Rho varying")
looic_cd4_list2 <- c(439.53, 461.99, 453.68, 423.97, 492.44, 447.03, 427.80)
looic_cd8_list2 <- c(492.32, 466.17, 493.95, 468.98, 475.66, 496.39, 482.38)

cd4_table2 <- model_compare(models_list2, looic_cd4_list2)
cd8_table2 <- model_compare(models_list2, looic_cd8_list2)

export_table_cd4 <- formattable(rbind(cd4_table2),
                             align = c('l', 'l', 'c', 'c'),
                             list(`Model` =
                                    formatter("span", style = ~ style(font.size='16px')),
                                  `Akaike weight %` =
                                    formatter("span", style = ~ style(font.size='16px')),
                                  `ΔLooIC` =
                                    formatter("span", style = x ~ ifelse(x == 0.00,
                                              style(color = customGreen, font.size='16px'),
                                              ifelse(x <= 6.0, style(color = customOrange, font.size='16px'), 
                                                     ifelse(x <100, style(font.size='16px'), NA)))))
)

export_table_cd8 <- formattable(rbind(cd8_table2),
                                align = c('l', 'l', 'c', 'c'),
                                list(`Model` =
                                       formatter("span", style = ~ style(font.size='16px')),
                                     `Akaike weight %` =
                                       formatter("span", style = ~ style(font.size='16px')),
                                     `ΔLooIC` =
                                       formatter("span", style = x ~ ifelse(x == 0.00,
                                                                            style(color = customGreen, font.size='16px'),
                                                                            ifelse(x <= 6.0, style(color = customOrange, font.size='16px'), 
                                                                                   ifelse(x <100, style(font.size='16px'), NA)))))
)

export_table_cd4
export_table_cd8

width_df = 400 * ncol(cd4_table2)
height_df = 150 * nrow(cd4_table2)
res_df = 40 * nrow(cd4_table2)

png(file.path("out_fit", "loo_table_cd4.png"), width=width_df,height=height_df, pointsize = 12, res = res_df, bg = NA)
dev.off()




