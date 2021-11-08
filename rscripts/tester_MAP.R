rm(list = ls()); gc()
library(rstan)
library(tidyverse)

expose_stan_functions(file.path('stan_models/Full_chimera/MAP_asm_deltavar_cd4.stan'))

x_r <- 100
x_i <- 61
params <- c(3.567291e+04, 1.094443e-02, 1.021495e-01, 1.485309e-01)
local_p <- c(0)

ymean <- math_reduce(params, local_p, x_r, x_i)
