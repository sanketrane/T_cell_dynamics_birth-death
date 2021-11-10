functions{
  // spline 1
  // Timecourse of thymic CD4 SP population -- changes with time
  real sp_numbers(real time) {
    real t0 = 1.0;
    real dpt0 = time - t0;     // days post t0
    real value; real fit1;
    // spline fitted separately to the counts of thymic SP8 cells
    // parameters estimated from spline fit to the timecourse of counts of source compartment -- SP CD8
    real theta0  = 9E4;    real theta_f = 68.0;  real n = 3.0;   real X  = 25.0;  real q = 4.25;
    //best fitting spline
    fit1 = theta0 + (theta_f * dpt0^n) * (1 - ((dpt0^q)/((X^q) + (dpt0^q))));

    if(time < t0){
      value = fit1;
    } else {
      value = fit1;
    }
    return value;
  }

  // spline 2
  // proportions of ki67 hi cells in source -- varies with time
  real eps_spline(real time){
    real value;
    // parameters estimated from spline fit to the timecourse of counts of source compartment -- SP CD8
    real eps_0 = 0.24510453; real eps_f = 0.01559996; real A = 14.83715328;
    real eps5 = exp(- eps_f * (5 + A)) + eps_0;    // the value of ki prop at day 5
    real fit;

    //best fitting spline
    if (time <= 5){
      fit = eps5; //* exp(-0.0.004255249 * (t-5));
    } else {
      fit  = exp(- eps_f * (time + A)) + eps_0;
    }
    return fit;
  }

  // initial age distribution of cells exisiting at t0
  real g_age(real age, real[] parms) {
    real t0 = 1.0;
    real N0 = parms[1];
    real value;

    // Flat, normalised age-distribution of initial cells
    if(age >= 0 && age <= t0) {
      value = N0/t0;
    } else {
      value = 0;
    }
    return value;
  }

  // Total influx into the naive T cell compartment from the thymus (cells/day)
  real theta_spline(real time, real[] parms){
    real N0 = parms[1];
    real t0 = 1.0;
    real psi;  real value;

    psi = g_age(0.0, parms)/sp_numbers(t0);
    value = psi * sp_numbers(time);
    return value;
  }

  // Ki67 distribution within the thymic influx -- varies with time
  real ki_dist_theta(real ki, real time, real[] parms){
    real k_bar = 1/exp(1.0);
    real value;
    real t0 = 1.0;
    real r_ki_theta = 0.9;     // thresold ki intensity within the thymic cohort below

    if(ki >= 0 && ki < k_bar){
      value = (1 - eps_spline(time))/k_bar;
    } else if(ki >= k_bar && ki < r_ki_theta) {
      value = 0.0;
    } else if (ki >= r_ki_theta && ki <= 1.0){
      value = (eps_spline(time)/(1 - r_ki_theta));
    } else {
      value = 0.0;
    }
    return value;
  }

  //  Ki67 distribution of cells exisiting in the periphery at t0
  real ki_dist_init(real ki){
    real value;
    real r_ki_init  = 3.0;         // parameter that shapes ki distribution within the init cohort

    if(ki >= 0.0 && ki <= 1.0){
      value = exp(ki * r_ki_init)/((exp(r_ki_init) - 1)/r_ki_init);
    }  else {
      value = 0.0;
    }
    return value;
  }

  // rate of cell division depending on cell age
  real rho_age(real age, real[] parms){
    real rho   = parms[3];
    real value  = rho;  // doesnt cahnge with cell age or host age
    return value;
  }

  // function that calculates intgral of net loss rate --  solved analytically to speed up the numerical integration
  real lambda_integ(real lo_lim, real up_lim, real[] parms){
    real del0  = parms[2];
    real rho   = parms[3];
    real r_del = parms[4];

    real value = ((del0/r_del) * (exp(-r_del * lo_lim) - exp(-r_del * up_lim))) - rho * (up_lim - lo_lim);
    return value;
  }

  // function that calculates intgral of net process rate --  solved analytically to speed up the numerical integration
  real alpha_integ(real lo_lim, real up_lim, real[] parms){
    real del0  = parms[2];
    real rho   = parms[3];
    real r_del = parms[4];

    real value = ((del0/r_del) * (exp(-r_del * lo_lim) - exp(-r_del * up_lim))) + rho * (up_lim - lo_lim);
    return value;
  }

  // Cell age distribution of the initial cohort
  real Asm_init_age(real age, real time, real[] parms) {
    real t0 = 1.0;

    real value = g_age(t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
    return value;
   }

   // Cell age distribution of the total theta cohort
   real Asm_theta_age(real age, real time, real[] parms) {

     real value =  theta_spline(time - age, parms) * exp(- lambda_integ(0, age, parms));
     return value;
    }

  // Cell age distribution of the whole population
  real Asm_total_age(real age, real time, real[] parms){
    real value;
    real t0 = 1.0;

    if(age < (time - t0)) {
      value =  theta_spline(time - age, parms) * exp(- lambda_integ(0, age, parms));
      } else {
        value = g_age(t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
    }

    return value; //Asm_init_age(age, time, parms) + Asm_theta_age(age, time,parms);
  }

  real[] asm_total_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = Asm_total_age(age, time, parms);

    return {value};
  }

  real solve_total_counts(real[] parms) {
    int x_i[0];
    real value;
    // integrate_ode_rk45(function, y0, t0, t, theta, x_r, x_i);
    // N0 is the y0 (init cond) -- also parameter[1] here
    // time = parms[5];

    value = integrate_ode_rk45(asm_total_ode, {0.0}, 0.0, rep_array(parms[5], 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  // Vectorised function for total pool size
  real[] N_total_time(data real[] time, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[5];
   params[1:4] = parms[1:4];

   for (i in 1:ndim){
     params[5] = time[i];
     y_solve[i] = solve_total_counts(params);
   }
   return y_solve;
  }

  // init cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_init_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 1.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * (time - t0)) ){
      value = g_age(t0, parms) * ki_dist_init(ki * exp(beta * (time - t0))) * exp(beta * (time - t0)) * exp(- alpha_integ(age - time + t0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_init_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 1.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_spline(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time - age, parms) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki'
  real[] U_total_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[5];
    real age = parms[6];
    real t0 = 1.0;
    real value;

    if (age < time - t0){
      value = U_theta_ki_age(ki, age, time, parms);
    } else {
      value = U_init_ki_age(ki, age, time, parms);
    }
    return {value};
  }

  // integral across ki values -- kbar to 1.0
  real U_total_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_total_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // the age distribution
  real U_total_age(real age, real[] parms){
    real params[6];
    real value;
    params[1:5] = parms[1:5];
    params[6] = age;

    return U_total_at(params);
  }

  // The integrand function for age distribution at time 'time'
  real[] U_total_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[5];
    real value = U_total_age(age, parms);

    return {value};
  }

  // integral across age values -- 0.0 to time
  real U_total_t(real[] parms){
    int x_i[0];

    real y_solve = integrate_ode_rk45(U_total_ode, {0.0}, 0.0, rep_array(parms[5], 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_total_time(data real[] time, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[5];
   params[1:4] = parms[1:4];

   for (i in 1:ndim){
     params[5] = time[i];
     y_solve[i] = U_total_t(params);
   }
   return y_solve;
  }

  vector lp_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_r)/2;
    int data_time[n] = x_i[1:n];                 // time -- host age
    real y_counts[n] = x_r[1:n];                 // counts data
    real y_ki[n] = x_r[(n+1):(2*n)];             // ki67 proportions data
    real y1_mean[n];                             // counts data predictions from the model
    real y2_mean[n];                             // ki67 proportions predictions from the model

    //PDE solution
    real counts_mean[n];
    real ki_counts[n];
    real kiprop_mean[n];

    //params
    real sigma_counts = global_params[5];
    real sigma_ki = global_params[6];

    //log likelihood
    real lp_manual;

    // unique timepoints in the data for ode solver and map of the indices of the unique timepoints
    int time_uniques[n];
    int time_indices[n];
    int already_exists;

    // first element is always unique
    int number_uniques = 1;
    time_uniques[1] = data_time[1];
    time_indices[1] = 1;

    for (i in 2:n) {
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

    // PDE solution -- predictions for total counts
    counts_mean[1:number_uniques] = N_total_time(to_array_1d(to_vector(time_uniques[1:number_uniques])/1.0), to_array_1d(global_params));
    ki_counts[1:number_uniques] = U_total_time(to_array_1d(to_vector(time_uniques[1:number_uniques])/1.0), to_array_1d(global_params));

    for (i in 1:number_uniques){
      kiprop_mean[i] = ki_counts[i]/counts_mean[i];
    }

    for (i in 1:n){
      y1_mean[i] = counts_mean[time_indices[i]];
      y2_mean[i] = kiprop_mean[time_indices[i]];
    }

    // joint log-likelihood since counts and ki67 are IID
    lp_manual = normal_lpdf(y_counts | to_vector(log(y1_mean)), sigma_counts) +
                normal_lpdf(y_ki | to_vector(logit(y2_mean)), sigma_ki);

    // return vector of joint likelihood
    return [lp_manual]';
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
  real<lower=1E4, upper=2E6> N0;
  real<lower=0.001, upper=0.5> delta;
  real<lower=0.001, upper=0.5> rho;
  real r_del;

  real<lower=0> sigma_counts;
  real<lower=0> sigma_ki;
}

transformed parameters{
  vector[6] global_params;
  global_params[1] = N0;
  global_params[2] = delta;
  global_params[3] = rho;
  global_params[4] = r_del;
  global_params[5] = sigma_counts;
  global_params[6] = sigma_ki;
}

model{
  N0 ~ normal(1E5, 3E4);
  delta ~ normal(0.05, 0.2);
  rho ~ normal(0.005, 0.2);
  r_del ~ normal(0.0, 0.2);

  sigma_counts ~ normal(0, 2);
  sigma_ki ~ normal(0, 2);

  target += sum(map_rect(lp_reduce, global_params, local_params, x_r, x_i));

  //log(counts) ~ normal(log(y1_mean), sigma_counts);
  //logit(ki_prop) ~ normal(logit(y2_mean), sigma_ki);
}

generated quantities{
  real y1_mean_pred[numPred];
  real ki_counts_pred[numPred];             // PDE prediction for ki+ counts
  real y2_mean_pred[numPred];                   // PDE prediction for ki prop
  real counts_pred[numPred];
  real kiprop_pred[numPred];

  // log likelihoods
  //vector[numObs] log_lik_counts;
  //vector[numObs] log_lik_ki;

  // PDE solution -- predictions for total counts
  y1_mean_pred = N_total_time(ts_pred, to_array_1d(global_params));
  ki_counts_pred = U_total_time(ts_pred, to_array_1d(global_params));

  for (i in 1: numPred){
   y2_mean_pred[i] = ki_counts_pred[i]/y1_mean_pred[i];
   counts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma_counts));
   kiprop_pred[i] = inv_logit(normal_rng(logit(y2_mean_pred[i]), sigma_ki));
  }

  // calculating log likelihoods
  //for (i in 1:numObs) {
  //  log_lik_counts[i] = normal_lpdf(log(counts[i]) | log(y1_mean[i]), sigma_counts);
  //  log_lik_ki[i]     = normal_lpdf(logit(ki_prop[i]) | logit(y2_mean[i]), sigma_ki);
  //}
}
