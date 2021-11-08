functions{
 //spline 1
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

  // Total influx into the naive T cell compartment from the thymus (cells/day)
  real theta_spline(real time, real psi){
    real value;

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

  real[] shm_chi(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
    real psi = parms[1];
    real delta = parms[2];
    real rho = parms[3];

    real dydt[4];
    real beta  = 1/3.5;            //rate of loss of ki67 -- mature naive cells

    // age of BMT in each recipient
    real ageAtBMT = parms[4];

    // ki hi donor
    dydt[1] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * eps_spline(time) + rho * (2 * y[2] +  y[1]) - (beta + delta) * y[1];
    // ki lo donor
    dydt[2] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * (1 - eps_spline(time)) + beta * y[1] - (rho  + delta) * y[2];

    // ki hi host
    dydt[3] = theta_spline(time, psi) * (1-Chi_spline(time - ageAtBMT)) * eps_spline(time) + rho * (2 * y[4] +  y[3]) - (beta + delta) * y[3];
    // ki lo host
    dydt[4] = theta_spline(time, psi) * (1-Chi_spline(time - ageAtBMT)) * (1 - eps_spline(time)) + beta * y[3] - (rho  + delta) * y[4];


    return dydt;
  }

  // solving for total counts of thymic and peripheral naive Tregs at time of BMT assuming the youngest animal as the t0
  // these counts form the initial conditions for other recipients N(0) = N_host(0)
  real[] solve_init(real ageAtBMT,
    real[] init_cond,                          // initial conditions at BMT in the youngest mouse
    real[] parms){

      real ta = 40;                            // age at BMT for the youngest animal
      real y_init[2, 4];
      real params_init[4];

      params_init[1:3] = parms[1:3];
      params_init[4] = ta;

      y_init[1] = init_cond;                 // init conditions at the earliest BMT (i.e. in younegst animal)
      y_init[2] = to_array_1d(integrate_ode_rk45(shm_chi, init_cond, ta, rep_array(ageAtBMT, 1), params_init, {0.0}, {0}));

      return y_init[2];
  }

  real[,] solve_ode_chi(real[] data_time,            // time point of observation
    real ageAtBMT,
    real[] init_cond,
    real[] parms){

      int num_obs = size(data_time);
      real y_solve[num_obs, 4];
      real params[4];

      real y0[4];
      real init_tb[4];                         // init conditions at the mean age of BMT for the group

      //solution for the initial conditions at the mean age of BMT for the group
      y0 = solve_init(ageAtBMT, init_cond, parms);

      // init conditions at the BMT
      init_tb[1] = 0;                                           //at tbmt - # donor is zero
      init_tb[2] = 0;                                           //at tbmt - # donor is zero
      init_tb[3] = y0[1] + y0[3];                               //at tbmt - all ki67Hi cells are host
      init_tb[4] = y0[2] + y0[4];                               //at tbmt - all ki67Lo cells are host


      params[1:3] = parms[1:3];
      params[4] = ageAtBMT;                                           // age at BMT

      y_solve = integrate_ode_rk45(shm_chi, init_tb, ageAtBMT, data_time, params, {0.0}, {0});

      return y_solve;
    }


    real[] shm_ont(real t, real[] y, real[] parms, real[] rdata,  int[] idata) {
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
      return to_array_1d(integrate_ode_rk45(shm_ont, init_cond, t0, rep_array(ts, 1), parms, {0.0}, {0}));
     }

    real[,] solve_ode_ont(real[] solve_time, real[] init_cond, real[] parms){
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
