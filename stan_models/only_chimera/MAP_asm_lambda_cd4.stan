functions{
  // proportions of ki67 hi cells in source -- varies with time
  real eps_spline(real time){
    real value;
    //parameters estimated from spline fit to the timecourse of ki67 proportions of source compartment -- SP CD4
    real eps_0 = 0.30895840; real eps_f = 0.03794025 ; real A = -16.65435251;
    real fit = exp(- eps_f * (time + A)) + eps_0;
    return fit;
  }

  // initial age distribution of cells exisiting at t0
  real g_age(real age, real[] parms) {
    real t0 = 40.0;
    real N0 = parms[1];
    real p0 = parms[4];
    real theta0;  real value;

    theta0 = N0 * p0/(exp(t0*p0)-1);

    // age-distribution of initial cells at t0
    if(age >= 0 && age <= t0) {
      value = theta0 * exp(p0 * age);
    } else {
      value = 0;
    }

    return value;
  }

  // Total influx into the naive T cell compartment from the thymus (cells/day)
  real theta_spline(real time, real[] parms){
    real N0 = parms[1];
    real p0 = parms[4];
    real t0 = 40.0;
    real theta0;  real value;

    //psi = g_age(0.0, parms)/sp_numbers(t0);
    theta0 = N0 * p0/(exp(t0*p0)-1);
    value = theta0 * exp(-0.003849627 * time);
    return value;
  }

  // influx of donor cells into the naive donor T cell compartment from the thymus (cells/day)
  real theta_donor(real time, real[] parms){
    real value;
    real tC = 42;
    real m = 2;

    //value = theta_spline(time, parms) * (time/tC)^m;
    if (time <= tC) {
      value = theta_spline(time, parms) * (time/tC)^m;
    } else {
      value = theta_spline(time, parms);
    }
    return value;
  }

  // influx of host derived cells into the naive host T cell compartment from the thymus (cells/day)
  real theta_host(real time, real[] parms){

    return theta_spline(time, parms) - theta_donor(time, parms);
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

  // function that calculates intgral of net process rate --  solved analytically to speed up the numerical integration
  real alpha_integ(real lo_lim, real up_lim, real[] parms){
    real lambda0 = parms[2];
    real r_lambda = parms[3];
    real rho   = parms[4];

    real value = lambda_integ(lo_lim, up_lim, parms) + 2 * rho * (up_lim - lo_lim);
    return value;
  }


  // Cell age distribution of the initial cohort
  real Asm_init_age(real age, real time, real[] parms) {
    real t0 = 40.0;

    real value = g_age(age-time+t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
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
    real t0 = 40.0;

    if(age < (time - t0)) {
      value =  theta_spline(time - age, parms) * exp(- lambda_integ(0, age, parms));
      } else {
        value = g_age(age-time+t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
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
  real[] N_total_time(real[] time, real[] parms){
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
  real[] N_pooled_time(real[] time, real[] tBMT, real[] parms){
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
  real[] N_host_time(real[] time, real[] tBMT, real[] parms){
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
  real[] N_donor_time(real[] time, real[] tBMT, real[] parms){
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
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * (time - t0)) ){
      value = g_age(age - time + t0, parms) * ki_dist_init(ki * exp(beta * (time - t0))) * exp(beta * (time - t0)) * exp(- alpha_integ(age - time + t0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_init_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_spline(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time - age) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_total_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
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
    real tBMT = parms[6];
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
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_host(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time - age) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Host_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // init cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real donor_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_donor(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time - age) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Donor_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }

    return value;
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki'
  real[] U_total_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real age = parms[7];
    real t0 = 40.0;
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
    real time = parms[7];
    real age = parms[8];
    real tBMT = parms[6];
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
    real time = parms[7];
    real age = parms[8];
    real tBMT = parms[6];
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
    real time = parms[7];
    real age = parms[8];
    real tBMT = parms[6];
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
    real params[7];
    real value;
    params[1:6] = parms[1:6];
    params[7] = age;

    return U_total_at(params);
  }

  // the age distribution
  real U_Pooled_age(real age, real[] parms){
    real params[8];
    real value;
    params[1:7] = parms[1:7];
    params[8] = age;

    return U_Pooled_at(params);
  }

  // the age distribution
  real U_host_age(real age, real[] parms){
    real params[8];
    real value;
    params[1:7] = parms[1:7];
    params[8] = age;

    return U_host_at(params);
  }

  // the age distribution
  real U_donor_age(real age, real[] parms){
    real params[8];
    real value;
    params[1:7] = parms[1:7];
    params[8] = age;

    return U_donor_at(params);
  }

  // The integrand function for age distribution at time 'time'
  real[] U_total_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real value = U_total_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_Pooled_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real value = U_Pooled_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_host_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real value = U_host_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_donor_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real value = U_donor_age(age, parms);

    return {value};
  }

  // integral across age values -- 0.0 to time
  real U_total_t(real[] parms){
    int x_i[0];
    real time = parms[6];

    real y_solve = integrate_ode_rk45(U_total_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_Pooled_t(real[] parms){
    int x_i[0];
    real time = parms[7];

    real y_solve = integrate_ode_rk45(U_Pooled_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_host_t(real[] parms){
    int x_i[0];
    real time = parms[7];

    real y_solve = integrate_ode_rk45(U_host_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_donor_t(real[] parms){
    int x_i[0];
    real time = parms[7];

    real y_solve = integrate_ode_rk45(U_donor_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_total_time(real[] time, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[6];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = time[i];
     y_solve[i] = U_total_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_Pooled_time(real[] time, data real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[7];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = tBMT[i];
     params[7] = time[i];
     y_solve[i] = U_Pooled_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_host_time(real[] time, data real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[7];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = tBMT[i];
     params[7] = time[i];
     y_solve[i] = U_host_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_donor_time(real[] time, data real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[7];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = tBMT[i];
     params[7] = time[i];
     y_solve[i] = U_donor_t(params);
   }
   return y_solve;
  }


  vector math_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_i); // n = 1
    int dat_t0 = x_i[1];                          // time zero -- for chimeras age at BM

    real chi_counts_mean[n];
    real host_counts_mean[n];
    real donor_counts_mean[n];

    vector[2*n] y_mean_stacked;

    // each shard has a single datpoint so its unique ****
    // PDE solution for chimera dataset -- x_r = data time and x_i = time at BMT
    chi_counts_mean = N_pooled_time(x_r,  to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    host_counts_mean = N_host_time(x_r, to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    donor_counts_mean = N_donor_time(x_r, to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));

    y_mean_stacked[1] = chi_counts_mean[1];
    y_mean_stacked[2] = donor_counts_mean[1]/chi_counts_mean[1];

    return y_mean_stacked;
  }

  // functions for transformation of fractions in (0,a), where a >=1
  real logit_inverse(real x){
     real ans;

     ans = exp(x)/(1+exp(x));

     return ans;
   }
}

data{
  int<lower = 0> numOnt;
  int<lower = 1> numChi;
  real<lower = 0> chi_counts[numChi];
  real<lower = 0> N_donor_fraction[numChi];
  int<lower  = 1> numPred;
  real<lower = 0> ts_pred_ont[numPred];
  real<lower = 0> ts_pred_chi1[numPred];
  real<lower = 0> ts_pred_chi2[numPred];
  real<lower = 0> ts_pred_chi3[numPred];
  real<lower = 0> tb_pred1[numPred];
  real<lower = 0> tb_pred2[numPred];
  real<lower = 0> tb_pred3[numPred];
  int n_shards;
  int<lower = 0> dat_t0[n_shards];       // nshards = numOnt + numChi
  int<lower = 0> dat_time[n_shards];     // nshards = numOnt + numChi
}

transformed data{
  int x_i[n_shards, 1];         // each shard gets a single data point
  real x_r[n_shards, 1];        // each shard gets a single data point

  // empty set of per shard params
  vector[0] local_params[n_shards];  // shard specific params --  useful for hierarchical modelling

  // data split into shards
  for (s in 1:n_shards){
   x_i[s, 1] = dat_t0[s];                       // age at BMT split
   x_r[s, 1] = dat_time[s];                     // time split
  }
}

parameters{
  real<lower=1E4, upper=2E6> N0;                  // total cells counts at t0
  real<lower=0.001, upper=0.5> lambda0;
  real r_lambda;
  real p0;
  real<lower=0> sigma_chi_counts;
  real<lower=0> sigma_Nfd;
}

transformed parameters{
  vector[4] global_params;
  vector[numChi] y3_mean;               // PDE prediction for counts from chimera data
  vector[numChi] y4_mean;               // PDE prediction for Nfd from chimera data
  vector[(2*numChi)] y_mean_stacked;        // compliled output across all nodes

  global_params[1] = N0;
  global_params[2] = lambda0;
  global_params[3] = r_lambda;
  global_params[4] = p0;

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:numChi){
    y3_mean[i] = y_mean_stacked[(4*numOnt)+ 2*i - 1];
    y4_mean[i] = y_mean_stacked[(4*numOnt)+ 2*i - 0];
  }
}

model{
  N0 ~ normal(5E5, 1.5E5);
  lambda0 ~ normal(0.05, 0.2);
  r_lambda ~ normal(0.0, 0.2);

  sigma_chi_counts ~ normal(0, 2);
  sigma_Nfd ~ normal(0, 2);

  log(chi_counts) ~ normal(log(y3_mean), sigma_chi_counts);
  logit(N_donor_fraction) ~ normal(logit(y4_mean), sigma_Nfd);
}

generated quantities{
  real y_chi_pred1[numPred, 2];
  real y_chi_pred2[numPred, 2];
  real y_chi_pred3[numPred, 2];

  real y3_mean_pred1[numPred];  real y4_mean_pred1[numPred];
  real y3_mean_pred2[numPred];  real y4_mean_pred2[numPred];
  real y3_mean_pred3[numPred];  real y4_mean_pred3[numPred];

  real chicounts_pred1[numPred]; real Nfd_pred1[numPred];
  real chicounts_pred2[numPred]; real Nfd_pred2[numPred];
  real chicounts_pred3[numPred]; real Nfd_pred3[numPred];

  real host_counts_pred1[numPred]; real host_counts_pred2[numPred];  real host_counts_pred3[numPred];
  real donor_counts_pred1[numPred]; real donor_counts_pred2[numPred];  real donor_counts_pred3[numPred];

  // log likelihoods
  vector[numChi] log_lik_chi_counts;
  vector[numChi] log_lik_Nfd;

  // PDE solution -- predictions for total counts, Nfd
  y3_mean_pred1 = N_pooled_time(ts_pred_chi1,  tb_pred1, to_array_1d(global_params));
  y3_mean_pred2 = N_pooled_time(ts_pred_chi2,  tb_pred2, to_array_1d(global_params));
  y3_mean_pred3 = N_pooled_time(ts_pred_chi3,  tb_pred3, to_array_1d(global_params));

  host_counts_pred1 = N_host_time(ts_pred_chi1, tb_pred1, to_array_1d(global_params));
  host_counts_pred2 = N_host_time(ts_pred_chi1, tb_pred2, to_array_1d(global_params));
  host_counts_pred3 = N_host_time(ts_pred_chi1, tb_pred3, to_array_1d(global_params));

  donor_counts_pred1 = N_donor_time(ts_pred_chi1, tb_pred1, to_array_1d(global_params));
  donor_counts_pred2 = N_donor_time(ts_pred_chi2, tb_pred2, to_array_1d(global_params));
  donor_counts_pred3 = N_donor_time(ts_pred_chi3, tb_pred3, to_array_1d(global_params));

  for (i in 1:numPred){
    y4_mean_pred1[i] = donor_counts_pred1[i]/y3_mean_pred1[i];
    y4_mean_pred2[i] = donor_counts_pred2[i]/y3_mean_pred2[i];
    y4_mean_pred3[i] = donor_counts_pred3[i]/y3_mean_pred3[i];

    chicounts_pred1[i] = exp(normal_rng(log(y3_mean_pred1[i]), sigma_chi_counts));
    chicounts_pred2[i] = exp(normal_rng(log(y3_mean_pred2[i]), sigma_chi_counts));
    chicounts_pred3[i] = exp(normal_rng(log(y3_mean_pred3[i]), sigma_chi_counts));

    Nfd_pred1[i] = logit_inverse(normal_rng(logit(y4_mean_pred1[i]), sigma_Nfd));
    Nfd_pred2[i] = logit_inverse(normal_rng(logit(y4_mean_pred2[i]), sigma_Nfd));
    Nfd_pred3[i] = logit_inverse(normal_rng(logit(y4_mean_pred3[i]), sigma_Nfd));
}

  // calculating log likelihoods
  for (i in 1:numChi) {
    log_lik_chi_counts[i] = normal_lpdf(log(chi_counts[i]) | log(y3_mean[i]), sigma_chi_counts);
    log_lik_Nfd[i]        = normal_lpdf(logit(N_donor_fraction[i]) | logit(y4_mean[i]), sigma_Nfd);
  }
}
