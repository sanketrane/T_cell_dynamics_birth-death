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

    real[] shm(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {

     real psi    = parms[1];
     real delta  = parms[2];
     real rho    = parms[3];

     real tb = idata[1];
     real Beta = 3.5;

     real dkdt[4];

     dkdt[1] = (psi * theta_spline(t, psi) * Chi_spline(t-tb) * 0.3) + rho * (2 * k[2] + k[1]) - ((1/Beta) + delta) * k[1];
     dkdt[2] = (psi * theta_spline(t, psi) * Chi_spline(t-tb) * (1-0.3)) + (1/Beta) * k[1] - (rho + delta) * k[2];

     dkdt[3] = (psi * theta_spline(t, psi) * (1-Chi_spline(t-tb)) * 0.3) + rho * (2 * k[4] + k[3]) - ((1/Beta) + delta) * k[3];
     dkdt[4] = (psi * theta_spline(t, psi) * (1-Chi_spline(t-tb)) * (1-0.3)) + (1/Beta) * k[3] - (rho + delta) * k[4];

     return dkdt;
   }

   real[] foreach_init(real tb, real ta, real[] init_cond, real[] parms, real[] rdata, int x_i){

     return to_array_1d(integrate_ode_rk45(shm, init_cond, ta, rep_array(tb, 1), parms, rdata, rep_array(x_i, 1)));
   }

   real[,] solve_init(real[] tb_time, real[] init_cond, real[] parms, real[] rdata, int x_i, int num_tb){
    real y_init[num_tb, 4];

    y_init[1] = init_cond;
    for (i in 2:num_tb){
      y_init[i] = foreach_init(tb_time[i], tb_time[1], init_cond, parms, rdata, x_i);
    }

    return y_init;
   }

   real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms, real[] rdata, int x_i) {

     return to_array_1d(integrate_ode_rk45(shm, init_cond, t0, rep_array(ts, 1), parms, rdata, rep_array(x_i, 1)));
    }

   real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time, int[] tb_index, int num_tb){
   real y_hat[num_index, 4];
   real y0[num_tb, 4];
   real init_tb[4];

   y0 = solve_init(tb_time, init_cond, parms, rdata, 45, num_tb);

   init_tb[1] = init_cond[1];                                           //at tbmt - donor ki67Hi subset size is zero
   init_tb[2] = init_cond[2];                                           //at tbmt - donor ki67Lo subset size is zero

   for (i in 1:num_index){
     init_tb[3] = y0[tb_index[i], 1] + y0[tb_index[i], 3];              //at tbmt - all ki67Hi cells would be host
     init_tb[4] = y0[tb_index[i], 2] + y0[tb_index[i], 4];              //at tbmt - all ki67Lo cells would be host
     y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
   }

   return y_hat;
   }

   real[,] solve_ode_pred(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time){
   real y_hat[num_index, 4];
   real y0[2, 4];
   real init_tb[4];

   y0 = solve_init(tb_time, init_cond, parms, rdata, 45, 2);

   init_tb[1] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero
   init_tb[2] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero

   init_tb[3] = y0[2, 1] + y0[2, 3];                                  //at tbmt - all ki67Hi cells would be host
   init_tb[4] = y0[2, 2] + y0[2, 4];                                  //at tbmt - all ki67Hi cells would be host

   y_hat[1] = init_tb;
   for (i in 2:num_index){
     y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
   }

   return y_hat;
   }
 }
