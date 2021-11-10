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
    real pop_dens = counts_func(time);

    real value = lambda / (1 + (pop_dens/count_bar)^3);
    return value;
  }

  real[] ode_func(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
    real dydt[1];

    // labelled naive cd8 T cells
    dydt[1] = - lambda_dens(time, parms) * y[1];

    return dydt;
  }

  real N_total(real time, real t0, real N0, real[] parms) {
    real value;
    // solves the ode for each timepoint from t0
    if (time == t0){
      value = 10^N0;
    } else {
      value = integrate_ode_rk45(ode_func, rep_array(10^N0, 1), t0, rep_array(time, 1), parms, {0.0}, {0})[1,1];
    }
    return value;
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

  real t0_group[J];          // group level info on day 0 of observations

  int numPred;            // size of predictor variable for prediction
  real tseq1[numPred];    // predictor variable for prediction of cohort1
  real tseq2[numPred];    // predictor variable for prediction of cohort1
  real tseq3[numPred];    // predictor variable for prediction of cohort1
  real tseq4[numPred];    // predictor variable for prediction of cohort1
  real tseq5[numPred];    // predictor variable for prediction of cohort1
  real ts_pred[numPred];    // predictor variable for predictions other quantities
  }

transformed data{
  vector[N] y1_counts;
  y1_counts = log(naive_labelled);                    // transforming cell counts for fitting
}

parameters{
  real<lower=0>  N0[K];         // total cells counts at t0
  real<lower=0> lambda;
  real<lower=0> count_bar;

  //hyper parameters
  real<lower=0> mu_N0;
  real<lower=0> sigma_N0;

  // variance of the likelihood
  real<lower=0> sigma_counts;
}

transformed parameters{
  real y_hat[N];
  real parms[2];

  // parameters without random effects
  parms[1] = lambda;
  parms[2] = count_bar;

  for (i in 1:N){
    // random effects in N0
    // different t0 for each group
    y_hat[i] = N_total(age_anim[i], t0_group[group_id[i]], N0[mouse_id[i]], parms);
  }
}

model{
  //priors
  lambda ~ normal(0.05, 0.2);
  count_bar ~ normal(5E6, 1E6);
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

  for (i in 1:numPred){
    // group 1
    y1_mean_pred[i] = N_total(tseq1[i], 14, mean(N0[1:12]), parms);

    // group 2
    y2_mean_pred[i] = N_total(tseq2[i], 28, mean(N0[13:29]), parms);

    // group 3
    y3_mean_pred[i] = N_total(tseq3[i], 42, mean(N0[30:41]), parms);

    // group 4
    y4_mean_pred[i] = N_total(tseq4[i], 70, mean(N0[42:53]), parms);

    // group 5
    y5_mean_pred[i] = N_total(tseq5[i], 189, mean(N0[54:66]), parms);
  }

   // calculating log likelihoods
   for (i in 1:N) {
     log_lik_counts[i] = normal_lpdf(y1_counts[i] | log(y_hat[i]), sigma_counts);
   }
}
