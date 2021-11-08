
Y1pred <- as.data.frame(model_fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_ont)%>%
  filter(timeseries >=5)

Cpred <- as.data.frame(model_fit, pars = "ontcounts_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_ont)%>%
  filter(timeseries >=5)

Y2pred <- as.data.frame(model_fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_ont)%>%
  filter(timeseries >=5)

kipred <- as.data.frame(model_fit, pars = "ontki_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_ont)%>%
  filter(timeseries >=5)


#### Age BMT bin1
Y3pred_bin1 <- as.data.frame(model_fit, pars = "y3_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

Chicounts_pred_bin1 <- as.data.frame(model_fit, pars = "chicounts_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

Y4pred_bin1 <- as.data.frame(model_fit, pars = "y4_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

Nfd_pred_bin1 <- as.data.frame(model_fit, pars = "Nfd_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

Y5pred_bin1 <- as.data.frame(model_fit, pars = "y5_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

donorki_pred_bin1 <- as.data.frame(model_fit, pars = "donorki_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

Y6pred_bin1 <- as.data.frame(model_fit, pars = "y6_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

hostki_pred_bin1 <- as.data.frame(model_fit, pars = "hostki_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi1)%>%
  filter(timeseries >=5)

#### Age BMT Bin2
Y3pred_bin2 <- as.data.frame(model_fit, pars = "y3_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

Chicounts_pred_bin2 <- as.data.frame(model_fit, pars = "chicounts_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

Y4pred_bin2 <- as.data.frame(model_fit, pars = "y4_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

Nfd_pred_bin2 <- as.data.frame(model_fit, pars = "Nfd_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

Y5pred_bin2 <- as.data.frame(model_fit, pars = "y5_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

donorki_pred_bin2 <- as.data.frame(model_fit, pars = "donorki_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

Y6pred_bin2 <- as.data.frame(model_fit, pars = "y6_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

hostki_pred_bin2 <- as.data.frame(model_fit, pars = "hostki_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi2)%>%
  filter(timeseries >=5)

#### Age BMT bin3
Y3pred_bin3 <- as.data.frame(model_fit, pars = "y3_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)

Chicounts_pred_bin3 <- as.data.frame(model_fit, pars = "chicounts_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)

Y4pred_bin3 <- as.data.frame(model_fit, pars = "y4_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)

Nfd_pred_bin3 <- as.data.frame(model_fit, pars = "Nfd_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)

Y5pred_bin3 <- as.data.frame(model_fit, pars = "y5_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)

donorki_pred_bin3 <- as.data.frame(model_fit, pars = "donorki_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)

Y6pred_bin3 <- as.data.frame(model_fit, pars = "y6_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)

hostki_pred_bin3 <- as.data.frame(model_fit, pars = "hostki_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred_chi3)%>%
  filter(timeseries >=5)


kidonor_df <- data.frame("ageBMT_group1" = Y5pred_bin1$median,
                         "ageBMT_group2" = Y5pred_bin2$median,
                         "ageBMT_group3" = Y5pred_bin3$median) %>%
  gather(key = ageBMTbin, value = donor_ki) 

kihost_df <- data.frame("ageBMT_group1" = Y6pred_bin1$median, 
                        "ageBMT_group2" = Y6pred_bin1$median, 
                        "ageBMT_group3" = Y6pred_bin1$median) %>%
  gather(key = ageBMTbin, value = host_ki)

kidonor_lb <- data.frame("ageBMT_group1" = Y5pred_bin1$lb,
                         "ageBMT_group2" = Y5pred_bin2$lb,
                         "ageBMT_group3" = Y5pred_bin3$lb) %>%
  gather(key = ageBMTbin, value = donor_ki) 

kihost_lb <- data.frame("ageBMT_group1" = Y6pred_bin1$lb, 
                        "ageBMT_group2" = Y6pred_bin1$lb, 
                        "ageBMT_group3" = Y6pred_bin1$lb) %>%
  gather(key = ageBMTbin, value = host_ki)

kidonor_ub <- data.frame("ageBMT_group1" = Y5pred_bin1$ub,
                         "ageBMT_group2" = Y5pred_bin2$ub,
                         "ageBMT_group3" = Y5pred_bin3$ub) %>%
  gather(key = ageBMTbin, value = donor_ki) 

kihost_ub <- data.frame("ageBMT_group1" = Y6pred_bin1$ub, 
                        "ageBMT_group2" = Y6pred_bin1$ub, 
                        "ageBMT_group3" = Y6pred_bin1$ub) %>%
  gather(key = ageBMTbin, value = host_ki)

ki_pred_lb <- data.frame("ageBMTbin"= kidonor_lb$ageBMTbin,
                         "Donor" = kidonor_lb$donor_ki,
                         "Host" = kihost_lb$host_ki) %>%
  gather(-ageBMTbin, key = Popln, value = "lb")

ki_pred_ub <- data.frame("ageBMTbin"= kidonor_ub$ageBMTbin,
                         "Donor" = kidonor_ub$donor_ki,
                         "Host" = kihost_ub$host_ki) %>%
  gather(-ageBMTbin, key = Popln, value = "ub")

ki_pred_bb <- data.frame(ki_pred_lb, "ub" = ki_pred_ub$ub) %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3,
                           ts_pred_chi1, ts_pred_chi2, ts_pred_chi3)) %>%
  arrange(ageBMTbin)


kipred_median <- data.frame(kidonor_df, kihost_df) %>%
  select(-"ageBMTbin.1") %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3)) %>%
  gather(-c(age.at.S1K, ageBMTbin), key = Popln, value = ki_prop) 

kidonor_bb <- data.frame("ageBMTbin" = kidonor_lb$ageBMTbin,
                         "lb" = kidonor_lb$donor_ki, "ub" = kidonor_ub$donor_ki) %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3))


