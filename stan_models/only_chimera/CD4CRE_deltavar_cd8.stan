functions{
  // initial age distribution at host age t0
  real g_age_dist(real age, real[] parms){
    real value;
    real N0 = parms[1];
    real day_stamp = parms[5];     // day of tamoxifen treatment for each group
    real t0 = parms[6];           // first tipme point of observation for each group

    // assuming that thymus continues to export stamped cells upto day 6 after the tamoxifen treatment
    // the number of stamoed cells exported daily is constant and therefore
    // the cell age distribution of stamped cells at t0  is uniform (more or less barring the loss dependent on cell age till t0)

    if (age < t0 - day_stamp - 1){
      value = 0.0;
    } else if (age <= t0 - day_stamp){
      value = (10^N0)/2;
    } else {
      value = 0.0;
    }

    return value;
  }

  // rate of cell division depending on cell age
  real lambda_age(real age, real[] parms){
    real delta = parms[2];
    real rho0  = parms[3];
    real r_rho = parms[4];

    real value  = delta - rho0 * exp(r_rho * age);

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

  // Total cell age distribution at a given time
  real[] G_age_time(real age, real[] y, real[] parms, real[] rdata, int[] idata){
    real t0 = parms[6];
    real time = parms[7];
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
    real day_stamp = parms[5];
    real time = parms[7];
    int x_i[0];

    real value = integrate_ode_rk45(G_age_time, {0.0}, time - day_stamp - 6.0, rep_array(time - day_stamp, 1), parms,  {0.0}, x_i)[1, 1];
    return value;
  }

  // Vectorised function for total pool size
  real N_total_time(data real time, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[7];
    params[1:6] = parms[1:6];

    for (i in 1:ndim){
      params[7] = time;
      y_solve[i] = N_total(params));
    }
    return y_solve;
  }
}
