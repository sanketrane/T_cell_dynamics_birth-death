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
      value = N0;
    } else {
      value = integrate_ode_rk45(ode_func, rep_array(N0, 1), t0, rep_array(time, 1), parms, {0.0}, {0})[1,1];
    }
    return value;
   }
}
