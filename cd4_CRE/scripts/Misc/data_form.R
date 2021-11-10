#Read datafile 
naive_df <- read.csv(file = "data/naive_labelled.csv") 

max(grp1_df$mouse_id)

grp1_df <- naive_df %>%
  filter(group_id == 5) %>%
  group_by(age_anim)%>%
  summarise(mean_grp = mean(naive_labelled),
            lb = quantile(naive_labelled, 0.11),
            ub = quantile(naive_labelled, 0.89)) %>%
  bind_cols(group_id = 1)

grp2_df <- naive_df %>%
  filter(group_id == 2) %>%
  group_by(age_anim)%>%
  summarise(mean_grp = mean(naive_labelled),
            lb = quantile(naive_labelled, 0.11),
            ub = quantile(naive_labelled, 0.89)) %>%
  bind_cols(group_id = 2)

grp3_df <- naive_df %>%
  filter(group_id == 3) %>%
  group_by(age_anim)%>%
  summarise(mean_grp = mean(naive_labelled),
            lb = quantile(naive_labelled, 0.11),
            ub = quantile(naive_labelled, 0.89)) %>%
  bind_cols(group_id = 3)

grp4_df <- naive_df %>%
  filter(group_id == 4) %>%
  group_by(age_anim)%>%
  summarise(mean_grp = mean(naive_labelled),
            lb = quantile(naive_labelled, 0.025),
            ub = quantile(naive_labelled, 0.975)) %>%
  bind_cols(group_id = 4)

grp5_df <- naive_df %>%
  filter(group_id == 5) %>%
  group_by(age_anim)%>%
  summarise(mean_grp = mean(naive_labelled),
            lb = quantile(naive_labelled, 0.11),
            ub = quantile(naive_labelled, 0.89)) %>%
  bind_cols(group_id = 5)

newdf <- rbind(grp1_df, grp2_df, grp3_df, grp4_df, grp5_df)

ggplot(data= newdf,  aes(x = age_anim, y = mean_grp))+
  geom_point(aes(col = as.factor(group_id)),alpha=0.8, size =2.5) + 
  #geom_line(aes(col = as.factor(group_id)), alpha=0.8, size =1) + 
  geom_errorbar(aes(ymin = lb, ymax=ub, col = as.factor(group_id))) +
  scale_color_discrete(name = NULL, guide = guide_legend(nrow = 5)) +guides(fill = F) +
  scale_y_continuous(limits = c(5e3, 2e6), trans="log10", breaks=c(1e5, 1e6, 1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
  labs(x = "Mouse age (days)", title = "Counts of timestamped naive CD8 T cells", y=NULL)  + 
  myTheme 
