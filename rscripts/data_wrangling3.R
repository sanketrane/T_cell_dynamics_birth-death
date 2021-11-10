library(readxl)
library(tidyverse)

#setwd("~/Desktop/Git_repositories/ontogeny_naive_Tcells")

ki67_chimeras <- read_excel("Datafiles/UPDATED_master_doc.xlsx", sheet = 9) %>%
  select("mouse.ID", "age.at.S1K", "SP.4nai", "LN.4nai", "SP.8nai", "LN.8nai") %>%
  na.omit() %>%
  mutate(ki_prop_ln4 = LN.4nai/100,
         ki_prop_sp4 = SP.4nai/100,
         ki_prop_ln8 = LN.8nai/100,
         ki_prop_sp8 = SP.8nai/100) %>%
  select("age.at.S1K", contains("ki_prop")) 

counts_chimeras <- read_excel("Datafiles/UPDATED_master_doc.xlsx", sheet = 4) %>%
  select("mouse.ID", "age.at.S1K", "SP.4nai", "LN.4nai", "SP.8nai", "LN.8nai") %>%
  na.omit() %>%
  mutate(counts_ln4 = LN.4nai,
         counts_sp4 = SP.4nai,
         counts_ln8 = LN.8nai,
         counts_sp8 = SP.8nai) %>%
  select("age.at.S1K", contains("counts")) 


counts_wt <- read_excel("Datafiles/UPDATED_master_doc.xlsx", sheet = 10) %>%
  filter(expt.location == "UCL") %>% 
  select("mouse.id","expt.location", "age.at.S1K", "SP.4nai", "LN.4nai", "SP.8nai", "LN.8nai") %>%
  na.omit() %>%
  mutate(counts_ln4 = LN.4nai,
         #counts_sp4 = SP.4nai,
         counts_ln8 = LN.8nai#,
         #counts_sp8 = SP.8nai
         ) %>%
  filter(age.at.S1K <=300) %>%
  select("mouse.id", "age.at.S1K", contains("counts")) 

ki67_wt <- read_excel("Datafiles/UPDATED_master_doc.xlsx", sheet = 11) %>%
  filter(expt.location == "UCL") %>%
  select("mouse.id", "expt.location", "age.at.S1K", "SP.4nai", "LN.4nai", "SP.8nai", "LN.8nai") %>%
  na.omit() %>%
  mutate(ki_prop_ln4 = LN.4nai/100,
        # ki_prop_sp4 = SP.4nai/100,
         ki_prop_ln8 = LN.8nai/100#,
         #ki_prop_sp8 = SP.8nai/100
         ) %>%
  filter(age.at.S1K <=300) %>%
  select("mouse.id", "age.at.S1K", contains("ki_prop")) 

donor_ki67 <- read_excel("Datafiles/UPDATED_master_doc.xlsx", sheet = 7) %>%
  select("mouse.ID", "age.at.S1K", "SP.4nai", "LN.4nai", "SP.8nai", "LN.8nai") %>%
  na.omit()  %>%
  mutate(ki_prop_ln4 = LN.4nai/100,
         ki_prop_sp4 = SP.4nai/100,
         ki_prop_ln8 = LN.8nai/100,
         ki_prop_sp8 = SP.8nai/100) %>%
  select("age.at.S1K", contains("ki_prop")) 

host_ki67 <- read_excel("Datafiles/UPDATED_master_doc.xlsx", sheet = 8) %>%
  select("mouse.ID", "age.at.S1K", "SP.4nai", "LN.4nai", "SP.8nai", "LN.8nai") %>%
  na.omit() %>%
  mutate(ki_prop_ln4 = LN.4nai/100,
         ki_prop_sp4 = SP.4nai/100,
         ki_prop_ln8 = LN.8nai/100,
         ki_prop_sp8 = SP.8nai/100) %>%
  select("age.at.S1K", contains("ki_prop")) 

data_counts <- read_excel("Datafiles/UPDATED_master_doc.xlsx", sheet = 4)%>%
  select("mouse.ID", "age.at.S1K", "SP.4nai", "LN.4nai", "SP.8nai", "LN.8nai") %>%
  na.omit()

ki67_ln_onto <- read_excel("Datafiles/ontogeny_data.xlsx", sheet = 3) %>%
  filter(organ == "LN") %>% select("mouse.id", "age.at.S1K.days", "CD4.nai", "CD8.nai")%>%
  mutate(ki_prop_ln4 = CD4.nai/100,
         ki_prop_ln8 = CD8.nai/100,
         age.at.S1K = age.at.S1K.days) %>% arrange(age.at.S1K)%>%
  select("mouse.id", "age.at.S1K", contains("ln"))

ki67_sp_onto <- read_excel("Datafiles/ontogeny_data.xlsx", sheet = 3) %>%
  filter(organ == "SP") %>% select("mouse.id", "age.at.S1K.days", "CD4.nai", "CD8.nai")%>%
  mutate(ki_prop_sp4 = CD4.nai/100,
         ki_prop_sp8 = CD8.nai/100,
         age.at.S1K = age.at.S1K.days) %>% arrange(age.at.S1K)


counts_ln_onto <- read_excel("Datafiles/ontogeny_data.xlsx", sheet = 2) %>%
  filter(organ == "LN") %>% select("mouse.id", "age.at.S1K.days", "CD4.nai", "CD8.nai")%>%
  mutate(counts_ln4 = CD4.nai,
         counts_ln8 = CD8.nai,
         age.at.S1K = age.at.S1K.days) %>% arrange(age.at.S1K) %>%
  select("mouse.id", "age.at.S1K", contains("ln"))

counts_sp_onto <- read_excel("Datafiles/ontogeny_data.xlsx", sheet = 2) %>%
  filter(organ == "SP") %>% select("mouse.id", "age.at.S1K.days", "CD4.nai", "CD8.nai")%>%
  mutate(counts_sp4 = CD4.nai,
         counts_sp8 = CD8.nai,
         age.at.S1K = age.at.S1K.days) %>% arrange(age.at.S1K)


ggplot() + 
  geom_point(data = host_ki67, aes(x = age.at.S1K, y = ki_prop_ln4), alpha =0.7, col = "#FF0000", size =2) +
  geom_point(data = donor_ki67, aes(x = age.at.S1K, y = ki_prop_ln4), col = "#00A08A", size =2) + scale_x_log10()

ggplot() + 
  geom_point(data = host_ki67, aes(x = age.at.S1K, y = ki_prop_ln8), alpha =0.7, col = "#FF0000", size =2) +
  geom_point(data = donor_ki67, aes(x = age.at.S1K, y = ki_prop_ln8), col = "#00A08A", size =2) + scale_x_log10()
  #geom_point(data = ki67_sp_onto, aes(x = age.at.S1K, y = ki_prop_sp4), col = "#5BBCD6", size =2) +
  geom_point(data = ki67_chimeras, aes(x = age.at.S1K, y = ki_prop_ln4), alpha =0.7,  col = "#F2AD00", size =2) +
  #geom_point(data = ki67_chimeras, aes(x = age.at.S1K, y = ki_prop_sp4), alpha =0.7,  col = "#00A08A", size =2) +
  scale_y_log10(limits = c(0.005, 1)) +theme_bw()

ggplot() + 
  geom_point(data = ki67_ln_onto, aes(x = age.at.S1K, y = ki_prop_ln8), alpha =0.7, col = "darkblue", size =2) +
  geom_point(data = ki67_sp_onto, aes(x = age.at.S1K, y = ki_prop_sp8), alpha =0.7, col = "darkred", size =2) +
  #geom_point(data = ki67_chimeras, aes(x = age.at.S1K, y = ki_prop_ln8), alpha =0.7,  col = "darkred", size =2) +
  #geom_point(data = ki67_chimeras, aes(x = age.at.S1K, y = ki_prop_sp8), alpha =0.7,  col = "darkgreen", size =2) +
  scale_y_log10(limits = c(0.005, 1))


ggplot() + 
  geom_point(data = counts_ln_onto, aes(x = age.at.S1K, y = counts_ln8), alpha = 0.7, col = "darkblue", size =2) +
  #geom_point(data = counts_sp_onto, aes(x = age.at.S1K, y = counts_sp8), alpha = 0.7, col = "darkred", size =2) +
  geom_point(data = counts_wt, aes(x = age.at.S1K, y = counts_ln8), alpha = 0.7,  col = "red", size =2) +
  #geom_point(data = counts_wt, aes(x = age.at.S1K, y = counts_sp8), alpha = 0.7,  col = "red", size =2) +
  scale_y_log10()


p1 <- ggplot() + 
  geom_point(data = counts_ln_onto, aes(x = age.at.S1K, y = counts_ln4), col = "#FF0000", size =3) +
  geom_point(data = counts_wt, aes(x = age.at.S1K, y = counts_ln4), alpha = 0.8,  col = "#5BBCD6", size =3) +
  scale_y_log10()+ xlim(0, 300)

p2 <- ggplot() + 
  geom_point(data = counts_ln_onto, aes(x = age.at.S1K, y = counts_ln8), col = "#FF0000", size =3) +
  geom_point(data = counts_wt, aes(x = age.at.S1K, y = counts_ln8), alpha = 0.8,  col = "#5BBCD6", size =3) +
  scale_y_log10()+ xlim(0, 300)


p3 <- ggplot() + 
  geom_point(data = ki67_ln_onto, aes(x = age.at.S1K, y = ki_prop_ln4), col = "#FF0000", size =3) +
  geom_point(data = ki67_wt, aes(x = age.at.S1K, y = ki_prop_ln4), alpha = 0.8,  col = "#5BBCD6", size =3) +
  scale_y_log10() + xlim(0, 300)

p4 <- ggplot() + 
  geom_point(data = ki67_ln_onto, aes(x = age.at.S1K, y = ki_prop_ln8),  col = "#FF0000", size =3) +
  geom_point(data = ki67_wt, aes(x = age.at.S1K, y = ki_prop_ln8), alpha = 0.8,  col = "#5BBCD6", size =3) +
  scale_y_log10() +xlim(0, 300)



cowplot::plot_grid(p1, p2, p3, p4, nrow = 2)

write.csv(counts_wt, "Datafiles/counts_wt.csv", row.names = F)
write.csv(ki67_wt, "Datafiles/ki67_wt.csv", row.names = F)






















