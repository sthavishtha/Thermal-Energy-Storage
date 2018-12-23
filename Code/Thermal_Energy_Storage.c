/*
---------------------------------------------------------------------------------------------------
Course Project of "Fundamentals of CFD Methods"
Objective of the project : Numerical simulation (CFD) of a thermocline thermal energy storage operating 
at high temperatures
Application of thermocline thermal energy storage : Concentrated solar power plants and process industries
Working principle behind thermal energy storage : Two steps - Charging and discharging
Charging - when the storage medium gets heated by the hot fluid flowing from top to bottom
Discharging - when the storage medium gets cooled by the fluid flowing from bottom to top

Behind the code:
 - Governing equations employ a two phase model, solving the energy equations
 - Finite Volume Method is adopted on a grid with uniform spacing
 - Both explicit and implicit time integration methods have been used

More details can be found in the Report 
---------------------------------------------------------------------------------------------------
 Author of the code : Sthavishtha Bhopalam Rajakumar, MSc in Mechanical Engineering
 Code updated on : 13-12-2017
 --------------------------------------------------------------------------------------------------
 Copyright Sthavishtha; All rights reserved
 --------------------------------------------------------------------------------------------------
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "function_list.h"

#define M_PI 3.14159265358979323846

int main()
{
        //initialization and declaration of some variables
        double time_cycle_start = 0.0, k_waveno1, k_waveno2, k_waveno;
	//k_wavnoe1 for the solid - Part 4 MMS 
	//k_wavenoe2 for the fluid - Part 4 MMS
	
	double exergy_cin = 0., exergy_cout = 0., exergy_din = 0., exergy_dout = 0., exergy_eff, t_inc_out, ex_error = 1.,
	  exergy_eff_new, q_c, q_d;

        int status, iteration = 0, iteration_mid = 0, user_value;
	double err = 1.0 ; //error to note the convergence
	time_current = 0.0;
        ovs_solid = 0, ovs_liquid = 0; //Enable them to 1 (separately) for respective MMS
	
	FILE* fileID;
	fileID = fopen("status_storage.csv","a");
	
	printf("Do you want to continue with the default values (Enter 0) or "
            "enter the values yourself (Enter a non-zero integer)\n");
	scanf("%d",&user_value);

	user_input(user_value); 
	
	//Wave numbers for fluid (1) and solid (2)
	k_waveno = 2*M_PI*7/height;

	if(ovs_liquid == 1 && ovs_solid == 1)
	  k_waveno2 = 2*k_waveno;

	else
	  k_waveno2 = k_waveno;

	k_waveno1 = k_waveno;

	delta_x = height/(double)(grid_no); //grid spacing

	h_vf = h_v/(rho_f*C_f*epsilon); //vol. heat transfer coefficient of fluid
	h_vs = h_v/(rho_s*C_s*(1 - epsilon)); //vol. heat transfer coefficient of solid

	printf("%f\t %f\n",h_vf,h_vs);

	//Dynamic memory implementation of the arrays
	t_fluid = (double*)malloc(grid_no*sizeof(double));
	t_solid = (double*)malloc(grid_no*sizeof(double));
	t_fluid_new = (double*)malloc(grid_no*sizeof(double)); //for the new iteration step
	t_fluid_iter = (double*)malloc(grid_no*sizeof(double)); //required for noting conv. (MMS)
	t_solid_new = (double*)malloc(grid_no*sizeof(double)); //for the new iteration step
	t_solid_iter = (double*)malloc(grid_no*sizeof(double)); //required for noting conv. (MMS)

	initialize(k_waveno1,k_waveno2); //initializing the initial temperatures
	
	visualization(iteration); 
		
	//	while(err > 1e-7) //convergence criterion for MMS convergence 
			while( ex_error > 1e-5) //convergence criterion for full scale simulation
	{

	  exergy_cin = 0., exergy_cout = 0., exergy_din = 0., exergy_dout = 0.;

		for(int i = 0;i < time_cycles; i++)
		{
			time_cycle_start = 0.0;
			iteration_mid = 0.0;

		while(time_cycle_start < time_per_cycle)
		{
			time_cycle_start = time_cycle_start + time_step;	
			time_current = time_current + time_step; //the actual time of the simulation
			status = status_storage(time_cycle_start);

			//iteration_mid++;
			
			//liq_solve(1,k_waveno);
			//solid_solve(status,k_waveno);
			coupled_solve_implicit(status,k_waveno1,k_waveno2);

			//equates new and old fluid and solid temperatures to be available for next iteration
			for(int l=0; l < grid_no; l++)
			  {
			     t_fluid[l] = t_fluid_new[l];
			     t_solid[l] = t_solid_new[l];
			  }

			if(iteration_mid == (int)(time_charge/time_step))
			  t_inc_out = t_fluid[grid_no-1] - t_fluid_bc_right;

			if(iteration_mid == (int)(time_per_cycle/time_step))
			  q_c = q_time();
			
			else if(iteration_mid == (int)((time_charge + time_idle)/time_step))
			  q_d = q_time();
				  
			if(status == 1)
			  {
			    exergy_cin += exergy_flux(0);
			    exergy_cout += exergy_flux(grid_no - 1);
			  }
			
			else if(status == 2)
			  {
			    exergy_din += exergy_flux(0);
			    exergy_dout += exergy_flux(grid_no - 1);
			  } 
			    
			if((iteration_mid)%(int)(10000.0/time_step) == 0)
			  {
				printf("Current Time (sec) = %lf\n",time_current);
				visualization(time_current);
			  }

			iteration_mid++;
			
		}
			
		}

		//Computes the error - to note the convergence for a specific set of iterations
		if(iteration%2000 == 0)
		  {
		    err = conv_error();
		    printf("Cycle number  = %d\t, Error = %6.8e\n",iteration,err);
		  }

		//t_fluid_iter is used for determining the convergence
		for(int a = 0;a < grid_no; a++)
		  {
		    t_fluid_iter[a] = t_fluid_new[a];
		    t_solid_iter[a] = t_solid_new[a];
		  }

		//Prints the temparatures after specific number of iterations
		if(iteration%10000 == 0 && iteration!=0)
		  visualization(iteration);

		iteration++;	

		printf("%lf \t %lf \t %lf \t %lf\n", exergy_dout, exergy_din, exergy_cin, exergy_cout);
	
	if (iteration == 1)
	  exergy_eff = fabs(exergy_dout - exergy_din)/(exergy_cin - exergy_cout);

	else 
	  {
	    exergy_eff_new = fabs(exergy_dout - exergy_din)/(exergy_cin - exergy_cout);
	    ex_error = fabs(exergy_eff_new - exergy_eff);
	    exergy_eff = exergy_eff_new;
	    printf("Intermediate values \n");
	    printf("Error = %lf\t Exergy new eff = %lf\n", ex_error, exergy_eff_new);
	    printf("Capacity factor = %lf\n", fabs(q_c - q_d)/q_max());
	    printf("Temperature increase at outflow = %lf\nf", t_inc_out);
			
	    }
	}

//		visualization(iteration);

	//printing the status after reaching convergence - for MMS
	for(int i = 0;i < time_cycles; i++)
		{
			time_cycle_start = 0.0;

			while(time_cycle_start <= time_per_cycle)
			  {
			    time_cycle_start = time_cycle_start + time_step;	
			    time_current = time_current + time_step;
			    status = status_storage(time_cycle_start);
			    fprintf(fileID,"%f\t %d\n",time_current,status);
			  }
		}

	printf("Final values of Cycle : Exergy Efficiency, Temperature increase at outflow "
	       "and Capacity factor are:\n");
	printf("Temperature increase at outflow = %lf\n", t_inc_out);
	printf("Cycle Exergy efficiency = %lf\n", exergy_eff);
	printf("Capacity factor = %lf\n", fabs(q_c - q_d)/q_max());
	

	
	// error_compute_ovs_solid(k_waveno2); //For MMS
	//error_compute_ovs_liquid(k_waveno1);
	
	fclose(fileID);
	
	//freeing the virtual memory used 
	free(t_fluid);
	free(t_fluid_new);
	free(t_solid);
	free(t_solid_new);
	free(t_fluid_iter);
	free(t_solid_iter);

	return 0;
		
}
