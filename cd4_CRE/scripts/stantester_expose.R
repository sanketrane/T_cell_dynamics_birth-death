#rm(list = ls()); gc()
library(rstan)
library(tidyverse)

## importing ASM stan model
ModelName <- 'APHA_N0mouse_Sigmoid2'
stan_file <- paste0(ModelName, '.stan')
expose_stan_functions(file.path('stan_models', stan_file))


# vector of parameters
parstan5_vec <- c(out_table$mean[72], out_table$mean[75],
                 out_table$mean[1], mean(out_table$mean[55:67]))
parstan2_vec <- c(out_table$mean[69], out_table$mean[75],
                  out_table$mean[1], mean(out_table$mean[14:30]))

#parstan_vec <- c(5.5, 0.08, 0.0007, 0.006)
rdata1_vec <- c(175, 189)
rdata2_vec <- c(7, 21)

## prediction points
tseq1 <- seq(14, 120, length.out = 100)
tseq2 <- seq(28, 120, length.out = 100)
tseq3 <- seq(42, 100, length.out = 100)
tseq4 <- seq(70, 160, length.out = 100)
tseq5 <- seq(189, 280, length.out = 100)

## stan solution
sim5 <- sapply(tseq5,  N_total_time, parms = c(parstan5_vec, rdata1_vec))
sim2 <- sapply(tseq2,  N_total_time, parms = c(parstan2_vec, rdata2_vec))

## labels for individual age cohorts
age_cohorts <- c(rep(5, 100), rep(2, 100))
## predictions data frame
naive_sim_df <- data.frame("timeseries" = c(tseq5, tseq2),
                           "group_id" = age_cohorts,
                           "median" = c(sim5, sim2))
#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv") %>%
  filter(group_id == c(2, 5))

ggplot()+
  geom_point(data = naive_sim_df, aes(x = timeseries, y = median, col = as.factor(group_id)),
            size =1.2) +
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(group_id)), size =3) +
  #geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled), col = "black", size =1.2) +
  scale_color_discrete(name = "Group ID", guide = guide_legend(nrow = 7)) +
  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +guides(col = F) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
  myTheme 

####### -------- plotting group and mouse wise lines ------ #########
## group 1
ntotal1_all <- matrix(nrow = 100, ncol = 12)
for (i in 2:13){
  parstan_group <- c(mean(out_table$mean[i]), out_table$mean[70], out_table$mean[1], 1, 14)
  ntotal1_all[, i-1] <- sapply(tseq1,  N_total_time, parms = parstan_group)
}

colnames(ntotal1_all) <- seq(1, 12)

ntotal_all1_df <- as.data.frame(ntotal1_all) %>%
  bind_cols('age_anim' = tseq1,
            'group_id' = rep("Cohort-1", 100)) %>%
  gather(-c(age_anim, group_id), key = 'mouse_id', value = 'naive_labelled')
  
# greoup 2
ntotal2_all <- matrix(nrow = 100, ncol = 17)
for (i in 14:30){
  parstan_group <- c(out_table$mean[i], out_table$mean[71], out_table$mean[1], 7, 28)
  ntotal2_all[, i-13] <- sapply(tseq2,  N_total_time, parms = parstan_group)
}

colnames(ntotal2_all) <- seq(14, 30)

ntotal_all2_df <- as.data.frame(ntotal2_all) %>%
  bind_cols('age_anim' = tseq2,
            'group_id' = rep("Cohort-2", 100)) %>%
  gather(-c(age_anim, group_id), key = 'mouse_id', value = 'naive_labelled')

# group 3
ntotal3_all <- matrix(nrow = 100, ncol = 12)
for (i in 31:42){
  parstan_group <- c(mean(out_table$mean[i]), out_table$mean[72], out_table$mean[1], 28, 42)
  ntotal3_all[, i-30] <- sapply(tseq3,  N_total_time, parms = parstan_group)
}

colnames(ntotal3_all) <- seq(31, 42)

ntotal_all3_df <- as.data.frame(ntotal3_all) %>%
  bind_cols('age_anim' = tseq3,
            'group_id' = rep("Cohort-3", 100)) %>%
  gather(-c(age_anim, group_id), key = 'mouse_id', value = 'naive_labelled')

# group 4
ntotal4_all <- matrix(nrow = 100, ncol = 12)
for (i in 43:54){
  parstan_group <- c(out_table$mean[i], out_table$mean[73], out_table$mean[1], 56, 70)
  ntotal4_all[, i-42] <- sapply(tseq4,  N_total_time, parms = parstan_group)
}

colnames(ntotal4_all) <- seq(43, 54)

ntotal_all4_df <- as.data.frame(ntotal4_all) %>%
  bind_cols('age_anim' = tseq4,
            'group_id' = rep("Cohort-4", 100)) %>%
  gather(-c(age_anim, group_id), key = 'mouse_id', value = 'naive_labelled')

## group 5
ntotal5_all <- matrix(nrow = 100, ncol = 13)
for (i in 55:67){
  parstan_group <- c(out_table$mean[i], out_table$mean[74], out_table$mean[1], 175, 189)
  ntotal5_all[, i-54] <- sapply(tseq5,  N_total_time, parms = parstan_group)
}

colnames(ntotal5_all) <- seq(55, 67)

ntotal_all5_df <- as.data.frame(ntotal5_all) %>%
  bind_cols('age_anim' = tseq5,
            'group_id' = rep("Cohort-5", 100)) %>%
  gather(-c(age_anim, group_id), key = 'mouse_id', value = 'naive_labelled')


ntotal_all_df <- rbind(ntotal_all1_df, ntotal_all2_df, ntotal_all3_df,
                       ntotal_all4_df, ntotal_all5_df)

naive_df$group_id <- as.factor(naive_df$group_id)
levels(naive_df$group_id) <- c('Cohort-1', "Cohort-2", "Cohort-3", "Cohort-4", "Cohort-5")
naive_sim_df$group_id <- as.factor(naive_sim_df$group_id)
levels(naive_sim_df$group_id) <- c('Cohort-1', "Cohort-2")

ggplot()+
  geom_point(data = ntotal_all_df, aes(x = age_anim, y = naive_labelled, col = as.factor(mouse_id)),
           alpha = 0.2, size =1.2) +
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(mouse_id)), size =3) +
  #geom_line(data = naive_sim_df, aes(x = age_anim, y = naive_labelled), col = "black", size =1.2) +
  scale_color_discrete(name = "Group ID", guide = guide_legend(nrow = 7)) +
  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +guides(col = F) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
  myTheme 

ggplot()+
  geom_point(data = ntotal_all_df, aes(x = age_anim, y = naive_labelled, col = as.factor(mouse_id)),
             alpha = 0.25, size =1.2) +
  geom_point(data= naive_df,  aes(x = age_anim, y = naive_labelled,
                                  col = as.factor(mouse_id)), size =3) +
  scale_color_discrete(name = "Group ID", guide = guide_legend(nrow = 7)) +
  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +guides(col = F) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
  myTheme  +
  facet_wrap(~ as.factor(group_id), scales = 'free_x') + guides(col = F) +
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_rect(fill = NA, colour = NA))


#ggsave(paste0(ModelName, '_', Population, 'simplot2.pdf'), width = 6, height = 4.5, units = "in")




### R sol
ntotal_R <- mcmapply(Asm_total_Time, tseq, MoreArgs = list(par_vec), mc.cores = num_cores)
kpos_R <- mcmapply(Kpos_total_Time, tseq, MoreArgs = list(par_vec), mc.cores = num_cores)

cd4_df <- data.frame('time' = tseq,
                     'counts' = ntotal_R,
                     'ki' = kpos_R)
write.csv(cd4_df, 'cd4_artf.csv', row.names = F)

ggplot()+
  geom_point(aes(x = tseq, y = ntotal_vec), size = 3, col = 2) + 
  geom_point(aes(x = tseq, y = ntotal_R), size = 1.2, col = 1) + 
  scale_y_log10()


ggplot()+
  geom_point(aes(x = tseq, y = kpos_vec/ntotal_vec), size = 3, col = 2) + 
  geom_point(aes(x = tseq, y = kpos_R/ntotal_R), size = 1.2, col = 1) + ylim(0,1) +
  scale_x_log10()



kseq <- seq(0, 1, 0.01)
plot(k_dist_init(kseq) ~ kseq)
kivec <- sapply(kseq, ki_dist_init, parms = parstan_vec)
points(kivec ~ kseq, col = 2)

## make artficial data
set.seed(12)
ts_pred <- round(c(runif(15, 5, 55), runif(25, 60, 300)), 0)


ntotal_df <- sapply(ts_pred, N_total_time, parms = parstan_vec)
init_df <- sapply(ts_pred, kpos_time, parms = parstan_vec)


err_tot <- rnorm(length(ts_pred),8e5, 1e6)
err_init <- rnorm(length(ts_pred), 2e4, 2e4)


artf_df <- data.frame('time_pred' = ts_pred,
                      'total_counts' = ntotal_df + err_tot,
                      'init_counts' = (init_df + err_init)) %>% arrange(time_pred) 


ggplot()+
  geom_point(data = artf_df, aes(x = time_pred, y = total_counts), size =1.2, col = 2) + 
  geom_line(aes(x = tseq, y = ntotal_vec), size =1) + 
  scale_y_log10()


ggplot()+
  geom_point(data = artf_df, aes(x = time_pred, y = ki_prop), size =1.2, col = 2) + 
  geom_line(aes(x = tseq, y = init_vec), size =1) 

write.csv(artf_df, 'artf_data.csv', row.names = F)
