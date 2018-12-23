#ifndef FUNCTION_LIST_H
#define FUNCTION_LIST_H

/* @BRIEF lists down all the functions to be used in the main program

   @FUNC

   liq_solve() :  Solves the discretized equations of the fluid discrete model at every time instant. Its parameter takes in the status - charging/discharging/idle
   conv_error() : To compute the error after a certain time steps (for convergence to be achieved or not) in the discrete model 
   liq_solve_mms() : Solves the discrete model using MMS - considering extra source terms 
   status_storage() : To determine the status of the storage medium at any time instant - charging/discharging/idle phase. Its parameter is the time_cycle_start 
   error_compute_ovs() : To compute the error in discrete model using MMS 
   user_input() : function to accept the user input of parameters   
   visualization() :  Prints the fluid, solid temperatures at any time step - for visualization by the user in a plotting package. A file at every cycle (until the convergence) is obtained.

*/

extern double diameter, height, u_f,diff_solid,diff_liq;
extern int grid_no, time_cycles, ovs_liquid, ovs_solid; //ovs - order verification study - separate for liquid and solid
extern double t_fluid_init,t_solid_init,t_fluid_bc_left, t_fluid_bc_right, delta_x ,time_step, time_charge, time_discharge, time_idle, time_current, time_per_cycle ;
extern double* t_fluid, *t_solid, *t_fluid_new, *t_solid_new, *t_fluid_iter ,*t_solid_iter;
extern double rho_s, rho_f, C_f, C_s, epsilon, h_vs, h_vf, h_v;

void liq_solve(int, double);
void solid_solve(int, double);
double conv_error();
//void liq_solve_mms();
int status_storage(double);
void error_compute_ovs_solid(double);
void error_compute_ovs_liquid(double);
void user_input(int);
void visualization(int);
void gauss_eliimination(double*, double*);
void coupled_solve_implicit(int, double, double);
void calc_param_sim();
void initialize(double, double);
double q_max();
double q_time(); 
double exergy_flux(int);
#endif
