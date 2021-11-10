functions{
  // initial age distribution at host age t0
  real g_age_dist(real age, real[] parms){
    real value;
    real N0 = parms[1];
    real day_stamp = parms[4];     // day of tamoxifen treatment for each group
    real t0 = parms[5];           // first tipme point of observation for each group

    // assuming that thymus continues to export stamped cells upto day 6 after the tamoxifen treatment
    // the number of stamoed cells exported daily is constant and therefore
    // the cell age distribution of stamped cells at t0  is uniform (more or less barring the loss dependent on cell age till t0)

    if (age < t0 - day_stamp - 6){
      value = 0.0;
    } else if (age <= t0 - day_stamp){
      value = (10^N0)/6;
    } else {
      value = 0.0;
    }

    return value;
  }

  // Rate of loss -- varies with cell age
  // integrand function
  real[] lambda_age(real age, real[] y, real[] parms, real[] rdata, int[] idata) {
    real lambda0 = parms[2];
    real r_lambda = parms[3];

    real value = lambda0 * exp(- r_lambda * age);
    return {value};
  }

  // function that calculates intgral of net loss rate --  solved analytically to speed up the numerical integration
  real lambda_integ(real lo_lim, real up_lim, real[] parms){
    real lambda0 = parms[2];
    real r_lambda = parms[3];

    real value = (lambda0/r_lambda) * (exp(-r_lambda * lo_lim) - exp(-r_lambda * up_lim));
    return value;
  }

  // Total cell age distribution at a given time
  real[] G_age_time(real age, real[] y, real[] parms, real[] rdata, int[] idata){
    real t0 = parms[5];
    real time = parms[6];
    real value;

    if (age > time - t0){
      value = g_age_dist(age - time + t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
      //exp(- integrate_ode_rk45(lambda_age, {0.0}, age - time + t0, rep_array(age, 1), parms, {0.0}, {0})[1, 1]);
    } else {
      value = 0.0;
    }

    return {value};
  }

  // Total counts of the stamped cells
  real N_total(real[] parms){
    real day_stamp = parms[4];
    real time = parms[6];
    int x_i[0];

    real value = integrate_ode_rk45(G_age_time, {0.0}, time - day_stamp - 6.0, rep_array(time - day_stamp, 1), parms,  {0.0}, x_i)[1, 1];
    return value;
  }

  // Vectorised function for total pool size
  real N_total_time(data real time, real[] parms){
    real params[6];
    params[1:5] = parms[1:5];
    params[6] = time;

    return N_total(params);
  }
}

data{
  int<lower=0> K;      // animal level assignments
  int<lower=0> J;      // cohort level assignments
  int<lower=0> N;      // number of observations

  vector[N] naive_labelled;   // array of observations
  real age_anim[N];       // predictor variable

  int<lower=1> mouse_id[N];   // array of animal assignments
  int<lower=1> group_id[N];   // array of group assignments

  real day_stamp[J];         // group level info on the day of labelling
  real t0_group[J];          // group level info on day 0 of observations

  int numPred;            // size of predictor variable for prediction
  real tseq1[numPred];    // predictor variable for prediction of cohort1
  real tseq2[numPred];    // predictor variable for prediction of cohort1
  real tseq3[numPred];    // predictor variable for prediction of cohort1
  real tseq4[numPred];    // predictor variable for prediction of cohort1
  real tseq5[numPred];    // predictor variable for prediction of cohort1
  }

transformed data{
  vector[N] y1_counts;
  y1_counts = log(naive_labelled);                    // transforming cell counts for fitting
}

parameters{
  real<lower=0>  N0[K];         // total cells counts at t0
  real<lower=0> lambda[J];      // vector of rate of loss for each animal
  real r_lambda;                // rate with which lambda changes with cell age

  //hyper parameters
  real<lower=0> mu_N0;
  real<lower=0> sigma_N0;
  real<lower=0> mu_lambda;
  real<lower=0> sigma_lambda;

  // variance of the likelihood
  real<lower=0> sigma_counts;
}

transformed parameters{
  real y_hat[N];
  real parms[5];

  // parameter with no random effects
  parms[3] = r_lambda;

  for (i in 1:N){
    parms[1] = N0[mouse_id[i]];              // random effects in N0
    parms[2] = lambda[group_id[i]];          // random effects in lambda0

    parms[4] = day_stamp[group_id[i]];       // different day of labelling for different groups
    parms[5] = t0_group[group_id[i]];        // different t0 for each group

    y_hat[i] = N_total_time(age_anim[i], parms);
  }
}

model{
  //priors
  r_lambda ~ normal(0.01, 0.1);
  sigma_counts ~ normal(0, 2);

  //hyper priors
  mu_N0 ~ normal(5.5, 1);
  sigma_N0 ~ normal(0.2, 0.5);
  mu_lambda ~ normal(0.05, 0.1);
  sigma_lambda ~ normal(0.005, 0.05);

  //model definition
  N0 ~ normal(mu_N0, sigma_N0);
  lambda ~ normal(mu_lambda, sigma_lambda);
  y1_counts ~ normal(log(y_hat), sigma_counts);
}


generated quantities{
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  real y3_mean_pred[numPred];
  real y4_mean_pred[numPred];
  real y5_mean_pred[numPred];
  // log likelihoods
   vector[N] log_lik_counts;

  // different sets of parameters for each group
  real parms1[5]; real parms2[5]; real parms3[5]; real parms4[5]; real parms5[5];

  for (i in 1:numPred){
    // group 1
    parms1[1] = mean(N0[1:12]);   parms1[2] = lambda[1];  parms1[3] = r_lambda;  parms1[4] = day_stamp[1];  parms1[5] = t0_group[1];
    y1_mean_pred[i] = N_total_time(tseq1[i], parms1);

    // group 2
    parms2[1] = mean(N0[13:29]);   parms2[2] = lambda[2];  parms2[3] = r_lambda;  parms2[4] = day_stamp[2];  parms2[5] = t0_group[2];
    y2_mean_pred[i] = N_total_time(tseq2[i], parms2);

    // group 3
    parms3[1] = mean(N0[30:41]);   parms3[2] = lambda[3];  parms3[3] = r_lambda;  parms3[4] = day_stamp[3];  parms3[5] = t0_group[3];
    y3_mean_pred[i] = N_total_time(tseq3[i], parms3);

    // group 4
    parms4[1] = mean(N0[42:53]);   parms4[2] = lambda[4];  parms4[3] = r_lambda;  parms4[4] = day_stamp[4];  parms4[5] = t0_group[4];
    y4_mean_pred[i] = N_total_time(tseq4[i], parms4);

    // group 5
    parms5[1] = mean(N0[54:66]);   parms5[2] = lambda[5];  parms5[3] = r_lambda;  parms5[4] = day_stamp[5];  parms5[5] = t0_group[5];
    y5_mean_pred[i] = N_total_time(tseq5[i], parms5);
  }

     // calculating log likelihoods
   for (i in 1:N) {
     log_lik_counts[i] = normal_lpdf(y1_counts[i] | log(y_hat[i]), sigma_counts);
   }
}
