functions{
  // initial age distribution at host age t0
  real g_age_dist(real age, real[] parms, real[] rdata){
    real value;
    real N0 = parms[1];
    real day_stamp = rdata[1];     // day of tamoxifen treatment for each group
    real t0 = rdata[2];           // first tipme point of observation for each group

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
    real t0 = rdata[2];
    real time = parms[4];
    real value;
    //print(age, time);

    if (age > time - t0){
      value = g_age_dist(age - time + t0, parms, rdata) * exp(- lambda_integ(age - time + t0, age, parms));
      //exp(- integrate_ode_rk45(lambda_age, {0.0}, age - time + t0, rep_array(age, 1), parms, {0.0}, {0})[1, 1]);
    } else {
      value = 0.0;
    }

    return {value};
  }

  // Total counts of the stamped cells
  real N_total(real[] parms, real[] rdata){
    real day_stamp = rdata[1];
    real time = parms[4];
    int x_i[0];
    real ub = parms[4] - rdata[1];

    real value = integrate_ode_rk45(G_age_time, {0.0}, time - day_stamp - 6, rep_array(time, 1), parms, rdata, x_i)[1, 1];
    return value;
  }

  // Vectorised function for total pool size
  real N_total_time(data real time, real[] parms, real[] rdata){
    real params[4];
    params[1:3] = parms[1:3];
    params[4] = time;

    return N_total(params, rdata);
  }
}

data{
  int<lower=0> K;      // number of animals
  int<lower=0> J;      // numbe of cohorts/groups
  int<lower=0> N;      // number of observations

  vector[N] naive_labelled;   // array of observations
  real age_anim[N];       // predictor variable

  int<lower=1, upper=K> mouse_id[N];   // array of animal assignments
  int<lower=1, upper=J> group_id[N];   // array of cohort assignments

  real day_stamp[J];         // group level info on the day of labelling
  real t0_group[J];          // group level info on day 0 of observations

  int numPred;            // size of predictor variable for prediction
  real ts_pred[numPred];  // predictor variable for prediction
  }

transformed data{
  vector[N] y1_counts;

  y1_counts = log(naive_labelled);                    // transforming cell counts for fitting
}

parameters{
  real<lower=0>  N0[J];     // total cells counts at t0
  real lambda[K];      // vector of rate of loss for each animal

  //real<lower=0> N0;
  //real lambda;
  real r_lambda;

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
  real parms[3];
  real rdata[2];

  //parms[1] = N0;
  //parms[2] = lambda;
  parms[3] = r_lambda;

  for (i in 1:N){
    parms[1] = N0[group_id[i]];              // group level variation
    parms[2] = lambda[mouse_id[i]];          // animal level variation

    rdata[1] = day_stamp[group_id[i]];       // different day of labelling for different groups
    rdata[2] = t0_group[group_id[i]];        // different t0 for each group

    y_hat[i] = N_total_time(age_anim[i], parms, rdata);
  }
}

model{
  //priors
  //N0 ~ normal(5.5, 1);
  //lambda ~ normal(0.05, 0.2);
  r_lambda ~ normal(0.01, 0.1);
  sigma_counts ~ normal(0, 2);

  //hyper priors
  mu_N0 ~ normal(5.5, 1);
  sigma_N0 ~ normal(0.2, 0.5);
  mu_lambda ~ normal(0.05, 0.1);
  sigma_lambda ~ normal(0.005, 0.05);

  //model definition
  //N0 ~ normal(mu_N0, sigma_N0);
  lambda ~ normal(mu_lambda, sigma_lambda);
  y1_counts ~ normal(log(y_hat), sigma_counts);
}


generated quantities{
  //real y1_mean_pred[numPred];

  //for (i in 1:numPred){
  //  //parms[1] = N0[mouse_id[i]];
  //  parms[1] = N0[mouse_id[i]];
  //  y_hat[i] = N_total_time(age_anim[i], parms);
  //}

  // log likelihoods
   vector[N] log_lik_counts;
   // calculating log likelihoods
   for (i in 1:N) {
     log_lik_counts[i] = normal_lpdf(y1_counts[i] | log(y_hat[i]), sigma_counts);
   }
}
