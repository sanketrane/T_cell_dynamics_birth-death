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

    if(ki >= 0 && ki < k_bar){
      value = (1 - eps_spline(time))/k_bar;
    } else if (ki >= k_bar && ki <= 1.0){
      value = (eps_spline(time)/(1 - k_bar));
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
    real value  = rho;  // doesnt change with cell age or host age
    return value;
  }

  // rate of cell loss depending on cell and host age
  real delta_age(real age, real time, real[] parms){
    real del0   = parms[2];
    real r_del = parms[4];
    real q_h   = parms[5];

    real del0_hostage = del0 + (4*del0/(1+(time/q_h)^5)); // delta of cell age = 0 varies with time.

    real value  = del0_hostage * exp(-r_del * age);  //  changes with cell age
    return value;
  }

  real lambda_integrand(real x,  real xc, real[] parms, real[] x_r,  int[] x_i) {
    real time = parms[6];
    real age = parms[7];

    real value = delta_age(x, x - age + time, parms) - rho_age(x, parms);
    return value;
  }

  real alpha_integrand(real x,  real xc, real[] parms, real[] x_r,  int[] x_i) {
    real time = parms[6];
    real age = parms[7];

    real value = delta_age(x, x - age + time, parms) + rho_age(x, parms);
    return value;
  }

  // Cell age distribution of the initial cohort
  real Asm_init_age(real age, real time, real[] parms) {
    real t0 = 1.0;
    real params[7];
    real value;

    params[1:5] = parms[1:5];
    params[6] = time;
    params[7] = age;

    value = g_age(t0, params) * exp(- integrate_1d(lambda_integrand, age - time + t0, age, params, {0.0}, {0}, 1E-12));
    return value;
   }

  real Asm_theta_age(real age, real time, real[] parms) {
    real t0 = 1.0;
    real params[7];
    real value;

    params[1:5] = parms[1:5];
    params[6] = time;
    params[7] = age;

    value =  theta_spline(time - age, params) * exp(- integrate_1d(lambda_integrand, 0, age, params, {0.0}, {0}, 1E-12));
    return value;
   }

 // Cell age distribution of the whole population
 real Asm_total_age(real age, real time, real[] parms){
   real value;
   real t0 = 1.0;

   if(age < (time - t0)) {
     value =  Asm_theta_age(age, time, parms);
     } else {
       value = Asm_init_age(age, time, parms);
   }

   return value; //Asm_init_age(age, time, parms) + Asm_theta_age(age, time,parms);
 }

 real[] asm_total_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
   real value;
   real time = parms[6];   // time (host age) input  as a param

   value = Asm_total_age(age, time, parms);

   return {value};
 }

 real solve_total_counts(real[] parms) {
   int x_i[0];
   real value;
   // integrate_ode_rk45(function, y0, t0, t, theta, x_r, x_i);
   // N0 is the y0 (init cond) -- also parameter[1] here

   value = integrate_ode_rk45(asm_total_ode, {0.0}, 0.0, rep_array(parms[6], 1), parms, {0.0}, x_i)[1, 1];
   return value;
 }

 // Vectorised function for total pool size
 real[] N_total_time(data real[] time, real[] parms){
  int ndim = size(time);
  real y_solve[ndim];
  real params[6];
  params[1:5] = parms[1:5];

  for (i in 1:ndim){
    params[6] = time[i];
    y_solve[i] = solve_total_counts(params);
  }
  return y_solve;
 }
   }
