functions{
 //spline 1
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

  real[] solve_unique(int[] data_time, real[] init_cond, real[] parms){
    int numObs = size(data_time);
    real y1_mean[numObs];                       // counts data predictions from the model
    real y2_mean[numObs];                       // ki67 data predictions from the model
    real y_hat[numObs, 2];                      // ode solution
    real y_mean_stacked[2 * numObs];

    // unique timepoints in the data for ode solver and map of the indices of the unique timepoints
    int time_uniques[numObs];
    int time_indices[numObs];
    int already_exists;

    // first element is always unique
    int number_uniques = 1;
    time_uniques[1] = data_time[1];
    time_indices[1] = 1;

    for (i in 2:numObs) {
      already_exists = 0;
      for (j in 1:number_uniques){
        if(time_uniques[j] == data_time[i]) {time_indices[i] = j; already_exists = 1;}
      }

      if(already_exists==0){
        number_uniques = number_uniques + 1;
        time_uniques[number_uniques] = data_time[i];
        time_indices[i] = number_uniques;
      }
    }

    // ODE solution -- predictions for the observed timecourse
    y_hat[1:number_uniques] = solve_ode(to_array_1d(to_vector(time_uniques[1:number_uniques])/1.0),
                              init_cond, parms);

    for (i in 1:numObs){
      y1_mean[i] = y_hat[time_indices[i], 1] + y_hat[time_indices[i], 2];
      y2_mean[i] = y_hat[time_indices[i], 1]/y1_mean[i];
    }

    for (i in 1:(numObs)){
        y_mean_stacked[i] = y1_mean[i];
        y_mean_stacked[i+numObs] = y2_mean[i];
      }

      return y_mean_stacked;
  }

  vector lp_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_r)/2;
    int data_time[n] = x_i[1:n];           // time -- host age
    real y_counts[n] = x_r[1:n];                 // counts data
    real y_ki[n] = x_r[(n+1):(2*n)];         // ki67 proportions data
    real y1_mean[n];                       // counts data predictions from the model
    real y2_mean[n];                       // ki67 data predictions from the model
    real y_mean_stacked[2 * n];                      // ode solution

    //params
    real sigma_counts = global_params[6];
    real sigma_ki = global_params[7];
    real N0 = global_params[4];
    real kappa0 = global_params[5];
    real init_cond[2];

    //log likelihood
    real lp_manual;

    // ODE solution -- predictions for the observed timecourse
    init_cond[1] = kappa0 * N0;
    init_cond[2] = (1 - kappa0) * N0;
    y_mean_stacked = solve_unique(data_time, init_cond, to_array_1d(global_params));

    // counts (y1_mean) and ki67 (y2_mean) proportions
    for (i in 1:n){
      y1_mean[i] = y_mean_stacked[i];
      y2_mean[i] = y_mean_stacked[i+n];
    }

    // joint log-likelihood since counts and ki67 are IID
    lp_manual = normal_lpdf(y_counts | to_vector(log(y1_mean)), sigma_counts) +
                normal_lpdf(y_ki | to_vector(logit(y2_mean)), sigma_ki);

    // return vector of joint likelihood
    return [lp_manual]';
  }

  vector map_tester(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_r)/2;
    int data_time[n] = x_i[1:n];           // time -- host age
    real y_counts[n] = x_r[1:n];                 // counts data
    real y_ki[n] = x_r[(n+1):(2*n)];         // ki67 proportions data
    real y1_mean[n];                       // counts data predictions from the model
    real y2_mean[n];                       // ki67 data predictions from the model
    real y_mean_stacked[2 * n];                      // ode solution

    //params
    real sigma_counts = global_params[6];
    real sigma_ki = global_params[7];
    real N0 = global_params[4];
    real kappa0 = global_params[5];
    real init_cond[2];

    //log likelihood
    real lp_manual;

    // ODE solution -- predictions for the observed timecourse
    init_cond[1] = kappa0 * N0;
    init_cond[2] = (1 - kappa0) * N0;
    y_mean_stacked = solve_unique(data_time, init_cond, to_array_1d(global_params));

    return to_vector(y_mean_stacked);
  }
}

data{
  int<lower = 1> numObs;
  int<lower = 0> data_time[numObs];
  real<lower = 0> counts[numObs];
  real<lower = 0> ki_prop[numObs];
  int<lower  = 1> numPred;
  real<lower = 0> ts_pred[numPred];
  int n_shards;
  }

transformed data{
  int M = numObs/n_shards;     // pershard numobs
  int x_i[n_shards, M];        // solve_time + time_index
  real x_r[n_shards, 2*M];     //2M because 2 variables, ki67 and counts

  // empty set of per shard params
  vector[0] local_params[n_shards];  // shard specific params --  useful for hierarchical modelling

  // data split into shards
  for (s in 1: n_shards){
    int i = 1 + (s-1) * M;                       // start index for ith shard
    int j = s * M;                               // end index for ith shard
    x_i[s, 1:M] = data_time[i:j];                // solve_time split
    x_r[s, 1:M] = log(counts[i:j]);              //counts split
    x_r[s, (M+1):(2*M)] = logit(ki_prop[i:j]);   //ki67 split
  }
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
  vector[7] global_params;
  vector[numObs] map_val;

  global_params[1] = psi;
  global_params[2] = delta;
  global_params[3] = rho;
  global_params[4] = N0;
  global_params[5] = kappa0;
  global_params[6] = sigma_counts;
  global_params[7] = sigma_ki;

  map_val = map_rect(map_tester, global_params, local_params, x_r, x_i);
  print(map_val);
}  

model{
  psi ~ normal(0.3, 0.2);
  delta ~ normal(0.03, 0.2);
  rho ~ normal(0.005, 0.2);
  N0 ~ normal(9E5, 1E5);
  kappa0 ~ normal(0.8, 0.1);

  sigma_counts ~ normal(0, 2);
  sigma_ki ~ normal(0, 2);

  target += sum(map_rect(lp_reduce, global_params,
                        local_params, x_r, x_i));
}

generated quantities{
  real y_hat_pred[numPred, 2];
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  real counts_pred[numPred];
  real kiprop_pred[numPred];
  real y_mean_stacked[2 * numObs];
  real y1_mean[numObs];
  real y2_mean[numObs];

  // log likelihoods
  vector[numObs] log_lik_counts;
  vector[numObs] log_lik_ki;
  vector[numObs] aic_ll;

  // initial conditions
  real init_cond[2];
  init_cond[1] = kappa0 * N0;
  init_cond[2] = (1 - kappa0) * N0;

  // ODE solution -- observed timecourse
  y_mean_stacked = solve_unique(data_time, init_cond, to_array_1d(global_params));

  for (i in 1:numObs){
    y1_mean[i] = y_mean_stacked[i];
    y2_mean[i] = y_mean_stacked[i+numObs];
  }

  // predictions for the whole timecourse
  y_hat_pred = solve_ode(ts_pred, init_cond, to_array_1d(global_params));

  for (i in 1:numPred){
    y1_mean_pred[i] = y_hat_pred[i, 1] + y_hat_pred[i, 2];
    y2_mean_pred[i] = y_hat_pred[i, 1] / y1_mean_pred[i];

    counts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma_counts));
    kiprop_pred[i] = inv_logit(normal_rng(logit(y2_mean_pred[i]), sigma_ki));
  }

  // calculating log likelihoods
  for (i in 1:numObs) {
    log_lik_counts[i] = normal_lpdf(log(counts[i]) | log(y1_mean[i]), sigma_counts);
    log_lik_ki[i]     = normal_lpdf(logit(ki_prop[i]) | logit(y2_mean[i]), sigma_ki);
    aic_ll[i]         = log_lik_counts[i] + log_lik_ki[i];
  }
}
