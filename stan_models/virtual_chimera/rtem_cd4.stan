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

  real Chi_spline( real time) {
    // chiEst is the level if stabilised chimerism in the source compartment
    // qEst is the rate with which cimerism chnages in the source compartment
    real chi;
    real chiEst = 0.847543332;
    real qEst = 0.050944623;

    if (time < 0){
      chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
    } else {
      chi = chiEst * (1 - exp(-qEst * time));
    }
    return chi;
  }

  real[] shm_chi(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
    real psi       = parms[1];
    real delta_nai = parms[2];
    real delta_rte = parms[3];
    real rho_nai   = parms[4];
    real rho_rte   = parms[5];
    real mu        = parms[6];

    real dydt[8];
    real beta  = 1/3.5;            //rate of loss of ki67 -- mature naive cells

    // age of BMT in each recipient
    real ageAtBMT = parms[7];

    // ki hi donor RTE
    dydt[1] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * eps_spline(time) + rho_rte * (2 * y[2] + y[1]) - (beta + delta_rte + mu) * y[1];
    // ki lo donor RTE
    dydt[2] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * (1 - eps_spline(time)) + beta * y[1] - (rho_rte + delta_rte + mu) * y[2];

    // ki hi donor mN
    dydt[3] = mu * y[1] + rho_nai * (2 * y[4] + y[3]) - (beta + delta_nai) * y[3];
    // ki lo donor mN
    dydt[4] = mu * y[2] + beta * y[3] - (rho_nai + delta_nai) * y[4];


    // ki hi host RTE
    dydt[5] = theta_spline(time, psi) * (1-Chi_spline(time - ageAtBMT)) * eps_spline(time) + rho_rte * (2 * y[6] + y[5]) - (beta + delta_rte + mu) * y[5];
    // ki lo host RTE
    dydt[6] = theta_spline(time, psi) * (1-Chi_spline(time - ageAtBMT)) * (1 - eps_spline(time)) +  beta * y[5] - (rho_rte + delta_rte + mu) * y[6];

    // ki hi mN
    dydt[7] = mu * y[5] + rho_nai * (2 * y[6] + y[7]) - (beta + delta_nai) * y[7];
    // ki lo mN
    dydt[8] = mu * y[6] + beta * y[7] - (rho_nai + delta_nai) * y[8];


    return dydt;
  }

  // solving for total counts of thymic and peripheral naive Tregs at time of BMT assuming the youngest animal as the t0
  // these counts form the initial conditions for other recipients N(0) = N_host(0)
  real[] solve_init(real ageAtBMT,
    real[] init_cond,                          // initial conditions at BMT in the youngest mouse
    real[] parms){

      real ta = 40;                            // age at BMT for the youngest animal
      real y_init[2, 8];
      real params_init[7];

      params_init[1:6] = parms[1:6];
      params_init[7] = ta;

      y_init[1] = init_cond;                 // init conditions at the earliest BMT (i.e. in younegst animal)
      y_init[2] = to_array_1d(integrate_ode_rk45(shm_chi, init_cond, ta, rep_array(ageAtBMT, 1), params_init, {0.0}, {0}));

      return y_init[2];
  }

  real[,] solve_ode_chi(real[] data_time,            // time point of observation
    real ageAtBMT,
    real[] init_cond,
    real[] parms){

      int num_obs = size(data_time);
      real y_solve[num_obs, 8];
      real params[7];

      real y0[8];
      real init_tb[8];                         // init conditions at the mean age of BMT for the group

      //solution for the initial conditions at the mean age of BMT for the group
      y0 = solve_init(ageAtBMT, init_cond, parms);

      // init conditions at the BMT
      init_tb[1] = 0;                                           //at tbmt - # donor is zero
      init_tb[2] = 0;                                           //at tbmt - # donor is zero
      init_tb[3] = 0;                                           //at tbmt - # donor is zero
      init_tb[4] = 0;                                           //at tbmt - # donor is zero
      init_tb[5] = y0[1] + y0[5];                               //at tbmt - all ki67Hi rte cells are host
      init_tb[6] = y0[2] + y0[6];                               //at tbmt - all ki67Lo rte cells are host
      init_tb[7] = y0[3] + y0[7];                               //at tbmt - all ki67Hi mN cells are host
      init_tb[8] = y0[4] + y0[8];                               //at tbmt - all ki67Lo mN cells are host

      params[1:6] = parms[1:6];
      params[7] = ageAtBMT;                                           // age at BMT

      y_solve = integrate_ode_rk45(shm_chi, init_tb, ageAtBMT, data_time, params, {0.0}, {0});

      return y_solve;
    }


    real[] shm_ont(real t, real[] y, real[] parms, real[] rdata,  int[] idata) {
      real psi       = parms[1];
      real delta_nai = parms[2];
      real delta_rte = parms[3];
      real rho_nai   = parms[4];
      real rho_rte   = parms[5];
      real mu        = parms[6];

      real dydt[4];
      real beta  = 1/3.5;            //rate of loss of ki67 -- mature naive cells

      // model that assumes that RTEs divide and die at different rates than mature naive T cells
      // ki hi RTE
      dydt[1] = theta_spline(t, psi) * eps_spline(t) + rho_rte * (2 * y[2] + y[1]) - (beta + delta_rte + mu) * y[1];
      // ki lo RTE
      dydt[2] = theta_spline(t, psi) * (1 - eps_spline(t)) + beta * y[1] - (rho_rte + delta_rte + mu) * y[2];
      // ki hi mN
      dydt[3] = mu * y[1] + rho_nai * (2 * y[4] + y[3]) - (beta + delta_nai) * y[3];
      // ki lo mN
      dydt[4] = mu * y[2] + beta * y[3] - (rho_nai + delta_nai) * y[4];

      return dydt;
    }

    real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms) {
      // solves the ode for each timepoint from t0
      return to_array_1d(integrate_ode_rk45(shm_ont, init_cond, t0, rep_array(ts, 1), parms, {0.0}, {0}));
     }

    real[,] solve_ode_ont(real[] solve_time, real[] init_cond, real[] parms){
      int num_solve = size(solve_time);
      real y_hat[num_solve, 4];
      // ode solution for the whole timecourse
      y_hat[1] = init_cond;
      for (i in 2:num_solve){
        y_hat[i] = foreach_ode(solve_time[i], solve_time[1], init_cond, parms);
      }
    return y_hat;
    }
   }
