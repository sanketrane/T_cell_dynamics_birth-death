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

  real Chi_spline( real time) {
    // chiEst is the level if stabilised chimerism in the source compartment
    // qEst is the rate with which cimerism chnages in the source compartment
    real chi;
    real chiEst = 0.85793486;
    real qEst = 0.04299229;

    if (time < 0){
      chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
    } else {
      chi = chiEst * (1 - exp(-qEst * time));
    }
    return chi;
  }

  // influx of donor cells into the naive donor T cell compartment from the thymus (cells/day)
  real theta_donor(real time, real[] parms){
    real value;
    real tBMT = parms[5];

    //value = theta_spline(time, parms) * (time/tC)^m;
    value = theta_spline(time, parms) * Chi_spline(time - tBMT);
    return value;
  }

  // influx of host derived cells into the naive host T cell compartment from the thymus (cells/day)
  real theta_host(real time, real[] parms){

    return theta_spline(time, parms) - theta_donor(time, parms);
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
    real r_rho = parms[4];

    real value  = rho * exp(-r_rho * age);

    return value;
  }

  // function that calculates intgral of net loss rate --  solved analytically to speed up the numerical integration
  real lambda_integ(real lo_lim, real up_lim, real[] parms){
    real delta = parms[2];
    real rho   = parms[3];
    real r_rho = parms[4];

    real value = (delta * (up_lim - lo_lim)) + ((rho/r_rho) * (exp(-r_rho * up_lim) - exp(-r_rho * lo_lim)));
    return value;
  }

  // function that calculates intgral of net process rate --  solved analytically to speed up the numerical integration
  real alpha_integ(real lo_lim, real up_lim, real[] parms){
    real delta = parms[2];
    real rho   = parms[3];
    real r_rho = parms[4];

    real value = (delta * (up_lim - lo_lim)) + ((rho/r_rho) * (exp(-r_rho * lo_lim) - exp(-r_rho * up_lim)));
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

    return value;
  }

  // Cell age distribution of the initial cohort
  real Asm_Host_init_age(real age, real time, real[] parms) {
    real tBMT = parms[5];
    real value;

    if (age >= time - tBMT){
      value = Asm_total_age(age - time + tBMT, tBMT, parms) * exp(- lambda_integ(age - time + tBMT, age, parms));
    } else {
      value = 0.0;
    }
    return value;
   }

   // Cell age distribution of the host theta cohort
   real Asm_Host_theta_age(real age, real time, real[] parms) {
     real tBMT = parms[5];
     real value;

     if (age < time - tBMT){
       value = theta_host(time - age, parms) * exp(- lambda_integ(0, age, parms));
     } else {
       value = 0.0;
     }
     return value;
   }

   // Cell age distribution of the donor theta cohort
   real Asm_Donor_theta_age(real age, real time, real[] parms) {
     real tBMT = parms[5];
     real value;

     if (age < time - tBMT){
       value = theta_donor(time - age, parms) * exp(- lambda_integ(0, age, parms));
     } else {
       value = 0.0;
     }
     return value;
   }

   // Cell age distribution of the total theta cohort
   real Asm_pooled_age(real age, real time, real[] parms) {
     real tBMT = parms[5];
     real value;

     if (age < time  - tBMT){
       value = Asm_Host_theta_age(age, time,  parms) + Asm_Donor_theta_age(age, time, parms);
     } else {
       value = Asm_Host_init_age(age, time, parms);
     }

     return value;
   }

   // Cell age distribution of the host cohort
   real Asm_host_age(real age, real time, real[] parms) {
     real tBMT = parms[5];
     real value;

     if (age < time  - tBMT){
       value = Asm_Host_theta_age(age, time,  parms);
     } else {
       value = Asm_Host_init_age(age, time, parms);
     }

     return value;
   }

   // Cell age distribution of the donor cohort
   real Asm_donor_age(real age, real time, real[] parms) {
     real tBMT = parms[5];
     real value;

     if (age < time - tBMT){
       value = Asm_Donor_theta_age(age, time, parms);
     } else {
       value = 0.0;
     }

     return value;
   }

  real[] Asm_total_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = Asm_total_age(age, time, parms);

    return {value};
  }

  real[] Asm_pooled_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[6];   // time (host age) input  as a param

    value = Asm_pooled_age(age, time, parms);

    return {value};
  }

  real[] Asm_host_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[6];   // time (host age) input  as a param

    value = Asm_host_age(age, time, parms);

    return {value};
  }

  real[] Asm_donor_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[6];   // time (host age) input  as a param

    value = Asm_donor_age(age, time, parms);

    return {value};
  }

  real solve_total_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[5];   // time (host age) input  as a param

    // integrate_ode_rk45(function, y0, t0, t, theta, x_r, x_i);
    value = integrate_ode_rk45(Asm_total_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_pooled_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[6];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_pooled_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_host_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[6];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_host_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_donor_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[6];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_donor_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
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

  // Vectorised function for total pool size
  real[] N_pooled_time(data real[] time, data real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[6];
    params[1:4] = parms[1:4];

    for (i in 1:ndim){
      params[5] = tBMT[i];
      params[6] = time[i];
      y_solve[i] = solve_pooled_counts(params);
    }
    return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_host_time(data real[] time, data real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[6];
    params[1:4] = parms[1:4];

    for (i in 1:ndim){
      params[5] = tBMT[i];
      params[6] = time[i];
      y_solve[i] = solve_host_counts(params);
    }
    return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_donor_time(data real[] time, data real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[6];
    params[1:4] = parms[1:4];

    for (i in 1:ndim){
      params[5] = tBMT[i];
      params[6] = time[i];
      y_solve[i] = solve_donor_counts(params);
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

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_total_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 1.0;
    real value;

    if (age < time - t0){
      value = U_theta_ki_age(ki, age, time, parms);
    } else {
      value = U_init_ki_age(ki, age, time, parms);
    }
    return value;
  }

  // init cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real host_init_ki_age(real ki, real age, real time, real[] parms){
    real tBMT = parms[5];
    real value;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions

    if (ki <= exp(-beta * (time - tBMT)) ){
      value = U_total_ki_age(ki * exp(beta * (time - tBMT)), age - time + tBMT, time, parms) * exp(beta * (time - tBMT)) * exp(- alpha_integ(age - time + tBMT, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Host_init_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real host_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 1.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_host(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time - age, parms) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Host_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // init cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real donor_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 1.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_donor(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time - age, parms) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Donor_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
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

  // The integrand function for age distribution of cells of ki intenisty 'ki' for the pooled donor and host compartments
  real[] U_Pooled_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real age = parms[7];
    real tBMT = parms[5];
    real value;

    if (age < time - tBMT){
      value = donor_theta_ki_age(ki, age, time, parms) + host_theta_ki_age(ki, age, time, parms);
    } else {
      value = host_init_ki_age(ki, age, time, parms);
    }
    return {value};
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki' for the pooled donor and host compartments
  real[] U_host_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real age = parms[7];
    real tBMT = parms[5];
    real value;

    if (age < time - tBMT){
      value = host_theta_ki_age(ki, age, time, parms);
    } else {
      value = host_init_ki_age(ki, age, time, parms);
    }
    return {value};
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki' for the pooled donor and host compartments
  real[] U_donor_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real age = parms[7];
    real tBMT = parms[5];
    real value;

    if (age < time - tBMT){
      value = donor_theta_ki_age(ki, age, time, parms);
    } else {
      value = 0.0;
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

  // integral across ki values -- kbar to 1.0 for the pooled donor and host compartments
  real U_Pooled_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_Pooled_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across ki values -- kbar to 1.0 for the pooled donor and host compartments
  real U_host_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_host_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across ki values -- kbar to 1.0 for the pooled donor and host compartments
  real U_donor_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_donor_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
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

  // the age distribution
  real U_Pooled_age(real age, real[] parms){
    real params[7];
    real value;
    params[1:6] = parms[1:6];
    params[7] = age;

    return U_Pooled_at(params);
  }

  // the age distribution
  real U_host_age(real age, real[] parms){
    real params[7];
    real value;
    params[1:6] = parms[1:6];
    params[7] = age;

    return U_host_at(params);
  }

  // the age distribution
  real U_donor_age(real age, real[] parms){
    real params[7];
    real value;
    params[1:6] = parms[1:6];
    params[7] = age;

    return U_donor_at(params);
  }

  // The integrand function for age distribution at time 'time'
  real[] U_total_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[5];
    real value = U_total_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_Pooled_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real value = U_Pooled_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_host_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real value = U_host_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_donor_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real value = U_donor_age(age, parms);

    return {value};
  }

  // integral across age values -- 0.0 to time
  real U_total_t(real[] parms){
    int x_i[0];
    real time = parms[5];

    real y_solve = integrate_ode_rk45(U_total_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_Pooled_t(real[] parms){
    int x_i[0];
    real time = parms[6];

    real y_solve = integrate_ode_rk45(U_Pooled_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_host_t(real[] parms){
    int x_i[0];
    real time = parms[6];

    real y_solve = integrate_ode_rk45(U_host_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_donor_t(real[] parms){
    int x_i[0];
    real time = parms[6];

    real y_solve = integrate_ode_rk45(U_donor_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
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

  // Vectorised function for theta div pool size
  real[] U_Pooled_time(data real[] time, data real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[6];
   params[1:4] = parms[1:4];

   for (i in 1:ndim){
     params[5] = tBMT[i];
     params[6] = time[i];
     y_solve[i] = U_Pooled_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_host_time(data real[] time, data real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[6];
   params[1:4] = parms[1:4];

   for (i in 1:ndim){
     params[5] = tBMT[i];
     params[6] = time[i];
     y_solve[i] = U_host_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_donor_time(data real[] time, data real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[6];
   params[1:4] = parms[1:4];

   for (i in 1:ndim){
     params[5] = tBMT[i];
     params[6] = time[i];
     y_solve[i] = U_donor_t(params);
   }
   return y_solve;
  }

  // functions for transformation of fractions in (0,a), where a >=1
  real[] asinsqrt_array(real[] x){
    int ndims = size(x);
    real answer[ndims];
    real a = 1.2;

    for (i in 1: ndims){
      answer[i] = asin(sqrt(x[i])/sqrt(a));
    }
    return answer;
  }

  real asinsqrt_real(real x){
    real a = 1.2;

    real answer = asin(sqrt(x)/sqrt(a));
    return answer;
   }

   real asinsqrt_inv(real x){
    real a = 1.2;

    real answer = a * (sin(x))^2;
    return answer;
   }
 }
