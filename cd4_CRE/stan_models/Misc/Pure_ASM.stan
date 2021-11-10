functions{
  // initial age distribution at host age t0
  real g_age_dist(real age, real[] parms, real[] rdata){
    real value;
    real N0 = parms[1];
    real day_stmp = rdata[1];     // day of tamoxifen treatment for each group
    real t0 = rdata[2];           // first tipme point of observation for each group

    // assuming that thymus continues to export stamped cells upto day 6 after the tamoxifen treatment
    // the number of stamoed cells exported daily is constant and therefore
    // the cell age distribution of stamped cells at t0  is uniform (more or less barring the loss dependent on cell age till t0)

    if (age < t0 - day_stmp - 6){
      value = 0.0;
    } else if (age <= t0 - day_stmp){
      value = N0/6;
    } else {
      value = 0.0;
    }

    return value;
  }

  // Rate of loss -- varies with cell age
  // integrand function
  real lambda_age(real age, real[] y, real[] parms, real[] rdata, int[] idata) {
    real lambda0 = parms[2];
    real r_lambda = parms[3];

    real value = lambda0 * exp(- r_lambda * age);
    return value;
  }

  // Total cell age distribution at a given time
  real G_age_time(real age, real[] y, real[] parms, real[] rdata, int[] idata){
    real t0 = rdata[2];
    real time = parms[4];

    real value = g_age_dist(age - time + t0, parms, rdata) *
    exp(- integrate_ode_rk45(lambda_age, {0.0}, age - time + t0, rep_array(age, 1), parms, rdata, {0})[1, 1]);

    return value;
  }

  // Total counts of the stamped cells
  real N_total(real[] parms, real[] rdata){
    real day_stmp = rdata[1];
    real time = parms[4];

    int x_i[0];

    real value = integrate_ode_rk45(G_age_time, {0.0}, time - day_stmp - 6, rep_array(time - day_stmp, 1), parms, rdata, x_i)[1, 1]
  }

  // Vectorised function for total pool size
  real[] N_total_time(data real[] time, real[] parms, real[] rdata){
   int ndim = size(time);
   real y_solve[ndim];
   real params[4];
   params[1:3] = parms[1:3];

   for (i in 1:ndim){
     params[4] = time[i];
     y_solve[i] = N_total(params);
   }
   return y_solve;
  }
}

data{
  int<lower=0> J;
  int<lower=0> numObs;
  int<lower=1, upper=J> mouse_id[numObs];
  vector[numObs] solve_time;
  vector[numObs] naive_labelled;
  vector[numObs] day_stmp;
  vector[numObs] t0_group;
  }

transformed data{
  vector[numObs] y1_counts;
  y1_counts = log(naive_labelled);                    // transforming cell counts for fitting

  real rdata[2];
  }

parameters{
  real N0[J];
  real lambda0[J];
  real r_lambda;

  real<lower=1E5, upper=2E6> mu_N0;                  // total cells counts at t0
  real mu_lambda0;

  real<lower=0> sigma_N0;
  real<lower=0> sigma_lamda0;
  real<lower=0> sigma_counts;
}

transformed parameters{
  vector[numObs] y_hat;

  real parms[3];

  for (j in 1:J){
    parms[1] = N0[j];
    parms[2] = lambda0[j];
    parms[3] = r_lambda;
  }


  for (i in 1:numObs){
    rdata[1] = day_stmp[i];
    rdata[2] = t0_group[i];

    y_hat[i] = N_total_time(solve_time, parms, rdata);
  }
}

model{
  for (j in 1:J){
    N0[j] ~ normal(mu_N0, sigma_N0);
    lambda0[j] ~ normal(mu_lambda0, sigma_lamda0);
  }
  mu_N0 ~ normal(1E6, 2.5E5);
  mu_lambda0 ~ normal(0.05, 0.2);
  r_lambda ~ normal(0.0, 0.2);

  sigma_counts ~ normal(0, 2);
  sigma_N0 ~ normal(0, 1e5);
  sigma_lamda0 ~ normal(0, 0.05);

  y1_counts ~ normal(log(y_hat), sigma_counts);
}
