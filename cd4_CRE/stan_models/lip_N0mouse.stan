functions{
  // spline for counts of the whole naive T cell population
  real counts_func(real time){
    real k = 7.06867;
    real n0 = 5.54431;
    real g = 0.15014;

    real value = 10^(k) * 10^(n0) * exp(g * time)/(10^(k) + 10^(n0) * (exp(g * time) - 1));
    return value;
  }

  // Rate of loss -- varies with cell density
  real lambda_dens(real time, real[] parms) {
    real lambda = parms[1];
    real count_bar = parms[2];

    // population density of the whole naive T cell compartment
    // calculated using the spline
    pop_dens = counts_func(time)

    real value = lambda * (1 + (pop_dens/count_bar)^3);
    return value;
  }

  // rate of loss of cells of age 0 varying with host age
  //real lambda_hostage(real time, real[] parms){
  //  real lambda0 = parms[2];
  //  real q_hostage = parms[4];
  //  real q2 = parms[5];

  //  return lambda0 * (1 + (q2/(1+(time/q_hostage)^3)));
  //}

  real[] ode(real t, real[] y, real[] parms, real[] rdata,  int[] idata) {
    real dydt[1];

    // labelled naive cd8 T cells
    dydt[1] = - lambda_dens(time, parms) * y[1];

    return dydt;
  }

  real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms) {
    // solves the ode for each timepoint from t0
    return to_array_1d(integrate_ode_rk45(ode, init_cond, t0, rep_array(ts, 1), parms, {0.0}, {0}));
   }

  real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms){
    int num_solve = size(solve_time);
    real y_hat[num_solve, 2];
    // ode solution for the whole timecourse
    y_hat[1] = init_cond;
    for (i in 2:num_solve){
      y_hat[i] = foreach_ode(solve_time[i], solve_time[1], init_cond, parms);
    }
  return y_hat;
 }

  // Total counts of the stamped cells
  real N_total(real[] parms){
    real day_stamp = parms[6];
    real time = parms[8];
    int x_i[0];

    real value = integrate_ode_rk45(G_age_time, {0.0}, time - day_stamp - 1.0, rep_array(time - day_stamp, 1), parms,  {0.0}, x_i)[1, 1];
    return value;
  }

  // Vectorised function for total pool size
  real N_total_time(data real time, real[] parms){
    real params[8];
    params[1:7] = parms[1:7];
    params[8] = time;

    return N_total(params);
  }
}

data{
  int<lower=0> K;      // animal level assignments
  int<lower=0> J;      // cohort level assignments
  int<lower=0> N;      // number of observations

  vector[N] naive_labelled;   // array of observations
  real age_anim[N];           // predictor variable

  int<lower=1> mouse_id[N];   // array of animal assignments
  int<lower=1> group_id[N];   // array of group assignments

  real day_stamp[J];         // group level info on the day of labelling
  real t0_group[J];          // group level info on day 0 of observations

  int numPred;            // size of predictor variable for prediction
  real tseq1[numPred];    // predictor variable for prediction of cohort1
  real tseq2[numPred];    // predictor variable for prediction of cohort2
  real tseq3[numPred];    // predictor variable for prediction of cohort3
  real tseq4[numPred];    // predictor variable for prediction of cohort4
  real tseq5[numPred];    // predictor variable for prediction of cohort5
  real ts_pred[numPred];    // predictor variable for predictions other quantities
  }

transformed data{
  vector[N] y1_counts;
  y1_counts = log(naive_labelled);                    // transforming cell counts for fitting
}

parameters{
  real<lower=0>  N0[K];         // total cells counts at t0
  real<lower=0> lambda;
  real<lower=0> r_lambda;
  real q_hostage;
  real q2;

  //hyper parameters
  real<lower=0> mu_N0;
  real<lower=0> sigma_N0;

  // variance of the likelihood
  real<lower=0> sigma_counts;
}

transformed parameters{
  real y_hat[N];
  real parms[7];

  // parameters without random effects
  parms[2] = 0.2631;
  parms[3] = r_lambda;
  parms[4] = 34;
  parms[5] = 1.4;

  for (i in 1:N){
    parms[1] = N0[mouse_id[i]];              // random effects in N0

    parms[6] = day_stamp[group_id[i]];       // different day of labelling for different groups
    parms[7] = t0_group[group_id[i]];        // different t0 for each group

    y_hat[i] = N_total_time(age_anim[i], parms);
  }
}

model{
  //priors
  //lambda ~ normal(0.027, 0.2);
  r_lambda ~ normal(0.01, 0.1);
  //q_hostage ~ normal(32, 2);
  //q2 ~ normal(1.3, 0.2);
  sigma_counts ~ normal(0, 2);

  //hyper priors
  mu_N0 ~ normal(5.5, 1);
  sigma_N0 ~ normal(0.2, 0.5);

  //model definition
  N0 ~ normal(mu_N0, sigma_N0);
  y1_counts ~ normal(log(y_hat), sigma_counts);
}


generated quantities{
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  real y3_mean_pred[numPred];
  real y4_mean_pred[numPred];
  real y5_mean_pred[numPred];
  real lambda0_hostage_pred[numPred];

  // log likelihoods
  vector[N] log_lik_counts;

  // different sets of parameters for each group
  real parms1[7]; real parms2[7]; real parms3[7]; real parms4[7]; real parms5[7];

  for (i in 1:numPred){
    // group 1
    parms1[1] = mean(N0[1:12]);   parms1[2] = lambda;  parms1[3] = r_lambda;   parms1[4] = q_hostage;  parms1[5] = q2;  parms1[6] = day_stamp[1];  parms1[7] = t0_group[1];
    y1_mean_pred[i] = N_total_time(tseq1[i], parms1);

    // group 2
    parms2[1] = mean(N0[13:29]);   parms2[2] = lambda;  parms2[3] = r_lambda;  parms2[4] = q_hostage;  parms2[5] = q2;  parms2[6] = day_stamp[2];  parms2[7] = t0_group[2];
    y2_mean_pred[i] = N_total_time(tseq2[i], parms2);

    // group 3
    parms3[1] = mean(N0[30:41]);   parms3[2] = lambda;  parms3[3] = r_lambda;  parms3[4] = q_hostage;  parms3[5] = q2;  parms3[6] = day_stamp[3];  parms3[7] = t0_group[3];
    y3_mean_pred[i] = N_total_time(tseq3[i], parms3);

    // group 4
    parms4[1] = mean(N0[42:53]);   parms4[2] = lambda;  parms4[3] = r_lambda;  parms4[4] = q_hostage;  parms4[5] = q2;  parms4[6] = day_stamp[4];  parms4[7] = t0_group[4];
    y4_mean_pred[i] = N_total_time(tseq4[i], parms4);

    // group 5
    parms5[1] = mean(N0[54:66]);   parms5[2] = lambda;  parms5[3] = r_lambda;  parms5[4] = q_hostage;  parms5[5] = q2;  parms5[6] = day_stamp[5];  parms5[7] = t0_group[5];
    y5_mean_pred[i] = N_total_time(tseq5[i], parms5);
  }

  //predictions for lambda0 varying with host age
  for (i in 1: numPred){
    lambda0_hostage_pred[i] = lambda_hostage(ts_pred[i], parms);
  }

   // calculating log likelihoods
   for (i in 1:N) {
     log_lik_counts[i] = normal_lpdf(y1_counts[i] | log(y_hat[i]), sigma_counts);
   }
}
