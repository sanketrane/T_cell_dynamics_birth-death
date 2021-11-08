functions{
  // spline 1
  // Timecourse of thymic CD4 SP population -- changes with time
  real sp_numbers(real time) {
    real t0 = 5.0;
    real dpt0 = time - t0;     // days post t0
    real value; real fit1;
    // spline fitted separately to the counts of thymic SP4 cells
    // parameters estimated from spline fit to the timecourse of counts of source compartment -- SP CD4
    real theta0  =  4.3E5;  real theta_f = 1.8E3;  real n = 2.1;   real X = 30.0;   real q = 3.7;
    //best fitting spline
    fit1 = theta0 + (theta_f * dpt0^n) * (1 - ((dpt0^q)/((X^q) + (dpt0^q))));

    if(time < t0){
      value = 0.0;
    } else {
      value = fit1;
    }
    return value;
  }

  // Total influx into the naive T cell compartment from the thymus (cells/day)
  real theta_spline(real time, real psi){
    real value;

    value = psi * sp_numbers(time);
    return value;
  }

  // spline2 --
   // proportions of ki67 hi cells in source -- varies with time
   real eps_spline(real time){
     real value;
     //parameters estimated from spline fit to the timecourse of ki67 proportions of source compartment -- SP CD4
     real eps_0 = 0.14965320; real eps_f = 0.03470231; real A = 3.43078629;
     //best fitting spline
     real fit = exp(-eps_f * (time + A)) + eps_0;

     return fit;
   }

  real[] shm(real t, real[] y, real[] parms, real[] rdata,  int[] idata) {
    real psi = parms[1];
    real delta = parms[2];
    real rho = parms[3];

    real dydt[2];
    real beta  = 1/3.5;            //rate of loss of ki67 -- mature naive cells

    // ki hi
    dydt[1] = theta_spline(t, psi) * eps_spline(t) + rho * (2 * y[2] +  y[1]) - (beta + delta) * y[1];
    // ki lo
    dydt[2] = theta_spline(t, psi) * (1 - eps_spline(t)) + beta * y[1] - (rho  + delta) * y[2];

    return dydt;
  }

  real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms) {
    // solves the ode for each timepoint from t0
    return to_array_1d(integrate_ode_rk45(shm, init_cond, t0, rep_array(ts, 1), parms, {0.0}, {0}));
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
  }

  data{
  int<lower  = 1> numObs;
  int<lower  = 1> num_index;
  int<lower  = 1> numPred;
  real<lower = 0> solve_time[num_index];
  real<lower = 0> ts_pred[numPred];
  real<lower = 0> counts[numObs];
  real<lower = 0> ki_prop[numObs];
  int<lower  = 0> time_index[numObs];
  }

  transformed data{
  real y_counts[numObs];
  real y_ki[numObs];
  real rdata[0];
  int idata[0];

  y_counts = log(counts);                 // transforming cell counts for fitting
  y_ki = logit(ki_prop);                  // transforming ki67 proportions for fitting
  }

  parameters{
    real<lower=1E4, upper=5E6> N0;                  // total cells counts at t0
    real<lower=0.001, upper=0.5> delta;
    real<lower=0.001, upper=0.5> rho;
    real<lower= 0, upper = 1> psi;
    real<lower= 0.2, upper = 1> kappa0;

    real<lower = 0> sigma_counts;
    real<lower = 0> sigma_ki;
  }

  transformed parameters{
  real y_hat[num_index, 2];
  real y1_mean[numObs];               // ODE prediction for counts
  real y2_mean[numObs];               // ODE prediction for ki prop

  real init_cond[2];
  real y1_0 = kappa0 * N0;
  real y2_0 = (1 - kappa0) * N0;

  real parms[3];
  parms[1] = psi;
  parms[2] = delta;
  parms[3] = rho;

  init_cond[1] = y1_0;
  init_cond[2] = y2_0;

  // ODE solution -- predictions for the observed timecourse
  y_hat = solve_ode(solve_time, init_cond, parms);

  for (i in 1:numObs){
    y1_mean[i] = y_hat[time_index[i], 1] + y_hat[time_index[i], 2];
    y2_mean[i] = y_hat[time_index[i], 1]/y1_mean[i];
  }
}

model{
  psi ~ normal(0.3, 0.2);
  delta ~ normal(0.03, 0.2);
  rho ~ normal(0.005, 0.2);
  N0 ~ normal(9E5, 1E5);
  kappa0 ~ normal(0.8, 0.1);

  sigma_counts ~ normal(0, 2);
  sigma_ki ~ normal(0, 2);

  y_counts ~ normal(log(y1_mean), sigma_counts);
  y_ki ~ normal(logit(y2_mean), sigma_ki);
}

generated quantities{
  real y_hat_pred[numPred, 2];
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  real counts_pred[numPred];
  real kiprop_pred[numPred];

  // log likelihoods
  vector[numObs] log_lik_counts;
  vector[numObs] log_lik_ki;
  vector[numObs] aic_ll;

  // ODE solution -- predictions for the whole timecourse
  y_hat_pred = solve_ode(ts_pred, init_cond, parms);

  for (i in 1:numPred){
    y1_mean_pred[i] = y_hat_pred[i, 1] + y_hat_pred[i, 2];
    y2_mean_pred[i] = y_hat_pred[i, 1] / y1_mean_pred[i];

    counts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma_counts));
    kiprop_pred[i] = inv_logit(normal_rng(logit(y2_mean_pred[i]), sigma_ki));
  }

  // calculating log likelihoods
  for (i in 1:numObs) {
    log_lik_counts[i] = normal_lpdf(y_counts[i] | log(y1_mean[i]), sigma_counts);
    log_lik_ki[i]     = normal_lpdf(y_ki[i] | logit(y2_mean[i]), sigma_ki);
    aic_ll[i]         = log_lik_counts[i] + log_lik_ki[i];
  }
}
