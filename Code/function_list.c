

#include "function_list.h"
#include<fstream>
#include<math.h>

#define M_PI 3.14159265358979323846

double diameter, height, u_f,diff_solid,diff_liq;
int grid_no, time_cycles, ovs_liquid, ovs_solid; //ovs - order verification study - separate for liquid and solid
double t_fluid_init,t_solid_init,t_fluid_bc_left, t_fluid_bc_right, delta_x,time_step, time_charge, time_discharge, time_idle, time_current, time_per_cycle ;
double* t_fluid, *t_solid, *t_fluid_new, *t_solid_new, *t_fluid_iter ,*t_solid_iter;
double rho_s, rho_f, C_f, C_s, epsilon, h_vs, h_vf, h_v, k_s, k_f, d_s, mass_flow, mu_f;


/**
   @brief Calculates the physical parameters for performing the full scale simulation
          Input to user_input() function
**/

void calc_param_sim()
{
  double h, h_fs, pr, re, nu_fs;

  //some extra parameters for Part 5 - for calc diff_solid/liq, u_f
  k_s = 2.0; //Thermal conductivity of the solid
  k_f = 0.52; //Thermal conductivity of the fluid
  mass_flow = 10; //mass flow rate of the fluid
  d_s = 0.03; //parameter used in calculation of overall heat transfer coefficient
  mu_f = 2.63; //fluid viscosity

  u_f = mass_flow*4/(rho_f*M_PI*epsilon*diameter*diameter); //fluid velocity

  pr = mu_f*C_f/k_f; //prandtl number of the fluid
  re = epsilon*rho_f*u_f*d_s/mu_f; //reynolds number 
  nu_fs = 0.255*pow(pr,0.333)*pow(re,0.667)/epsilon; //Nusselt number 
  h_fs = nu_fs*k_f/d_s; //fluid-solid heat transfer coefficient
  h = 1.0/((1.0/h_fs) + (d_s/(10.0*k_s))); //overall heat transfer coefficient
  h_v = 6.0*(1.0 - epsilon)*h/d_s; //volumetric heat transfer coefficient

  diff_solid = k_f/(epsilon*rho_f*C_f); //diffusivity of the solid
  diff_liq = k_s/((1-epsilon)*(rho_s*C_s)); //diffusivity of the fluid

}

/**
   @brief Accepts the user input data or uses the standard specified values edited by the programmer
   @param[in] Integer b to determine if the default or user specified values have to be accepted
**/

void user_input(int b)
{	
	
	if(b == 0) //standard specified values 
	{
	        //Geometrical parameters
	        diameter = 4; //diameter of the storage medium
		height = 23.8732; //height of the storage medium
		grid_no = 200; //number of grid cells

		//Temperatures - initial and boundary conditions
		t_fluid_init = 293; //initial fluid temperature
		t_solid_init = 293; //initial solid temperature
		t_fluid_bc_left = 873; //charging phase
		t_fluid_bc_right = 293; //discharging phase

		//fluid, solid dependent variables

		h_v = 1000;
		rho_s = 2600;  //density of the solid phase 
		rho_f = 1835.6; //density of the fluid phase
		C_s = 900.0; //specific heat capacity of the solid
		C_f = 1511.8; //specific heat capacity of the fluid

		u_f = 1e-5; 
		diff_solid = 9e-7; 
		diff_liq = 2e-7; 
		epsilon = 0.4; //factor to account for the mass flow (due to the
		//obstruction by the solid
  
               //Calculates the parameters for the simulation - uncommment this for MMS
	       // Used only for Part 5 - Validation
	        calc_param_sim();
				
		//cycle variables - time
		time_charge = 21600; //chanrging time
		time_discharge = 21600; //discharging time
		time_idle = 21600; //Idle time 
		time_cycles = 1; //Number of cycles
		time_step = 1e-2; //Delta t 
	}
				
	else
	{	
		printf("Numerical Simulation of a thermocline thermal Energy Storage - A CFD study\n");
		printf("Assuming the storage medium to be cylindrical, please input the diameter of the storage medium\n");
		scanf("%lf",&diameter);
		printf("Please input the height of the storage medium\n");
		scanf("%lf",&height);
		printf("Please input the number of cells in the grid - one dimensional grid\n");
		scanf("%d",&grid_no);
		printf("Please input the initial temperatures of the fluid or solid phases (which are initially taken to be equal)\n");
		scanf("%lf",&t_fluid_init);
		t_solid_init=t_fluid_init;
		printf("Please input the left temperature boundary condition for the fluid\n");
		scanf("%lf",&t_fluid_bc_left);
		printf("Please input the right temperature boundary condition for the fluid\n");
		scanf("%lf",&t_fluid_bc_right);
		printf("Please input the fluid velocity u_f\n");
		scanf("%lf",&u_f);
		printf("Please input the diffusivity of liquid\n");
		scanf("%lf",&diff_liq);
		printf("Please input the diffusivity of solid\n");
		scanf("%lf",&diff_solid);
		printf("For determination of the state of the storage medium\n ");
		printf("Enter the Charging time\n");
		scanf("%lf",&time_charge);
		printf("Enter the Discharging time\n");
		scanf("%lf",&time_discharge);
		printf("Enter the idle time\n");
		scanf("%lf",&time_idle);
		printf("Enter the total number of cycles\n");
		scanf("%d",&time_cycles);
	}
	
	time_per_cycle =  time_charge + 2*time_idle + time_discharge; //1 cycle time = sum of durations of all 4 time intervals
	
}

/**
   @brief Gaussian Elimination to solve the 2 equations - fluid and solid for point implicit method
   @param[in] Fluid temperature at one particular point
   @param[in] Solid temperature at one particular  point
**/

void gauss_elimination(double* fluid_val, double* solid_val)
{
	  double a,b,c,d;
	  a = 1 + (h_vf*time_step);
	  b = -h_vf*time_step;
	  c = -h_vs*time_step;
	  d = 1 + (h_vs*time_step);
	  *solid_val = (*solid_val - (c/a)*(*fluid_val))/(d - (c*b/a));
	  *fluid_val = (*fluid_val - b*(*solid_val))/a;
}


/**
   @brief Prints the solid and fluid temperature data into the file
   @param[in] Iteration number
**/

void visualization(int iter)
{
        FILE* fileID;
	char file_name[80];
	sprintf(file_name,"temperatures_%05d.dat",iter); 
	fileID = fopen(file_name,"w");
	fprintf(fileID,"Position\t Fluid Temperature\t Solid temperature\n");

	for(int i=0;i<grid_no;i++)
		fprintf(fileID,"%lf\t %lf\t %lf\n",0.5*(double)(2*i + 1)*delta_x,t_fluid[i],t_solid[i]);
	fclose(fileID);
}

/**
   @brief initializes the initial temperature values of fluid and solid
**/

void initialize(double k_waveno1, double k_waveno2)
{
  double angle;
  if (ovs_solid == 0) //if MMS for solid has been disabled
	  {
	    //Initialization for a normal simulation
	  for(int i = 0;i < grid_no; i++)
	    {
			t_solid[i] = t_solid_init;
			t_solid_new[i] = t_solid_init;
			t_solid_iter[i] = t_solid_init;
	    }
	  }
	else
	{ 
	  //Initialization for MMS
	  for(int i = 0;i < grid_no; i++)
	    {
	                angle = k_waveno2*delta_x*0.5*(double)(2*i + 1);
			
			//Initial conditions not exactly equal to angle -reason mentioned in the report
			t_solid[i] = 0.8*cos(angle);					
			t_solid_iter[i] = 0.8*cos(angle);
			t_solid_new[i] = 0.8*cos(angle);
	    }
	}


	if (ovs_liquid == 0)
	  {
	    //Initialization for a normal simulation	    
	    for(int i = 0;i < grid_no; i++)
	    {
	      t_fluid[i] = t_fluid_init;
	      t_fluid_new[i] = t_fluid_init;
	      t_fluid_iter[i] = t_fluid_init;
	     }
	  }
	else
	{
   	  //Initialization for MMS
	      for(int i = 0;i < grid_no; i++)
		{
		        angle = k_waveno1*delta_x*0.5*(double)(2*i + 1);

			//Initial conditions not exactly equal to angle -reason mentioned in the report
			t_fluid[i] = 0.8*cos(angle);
			t_fluid_iter[i] = 0.8*cos(angle);
			t_fluid_new[i] = 0.8*cos(angle);
		}
	}
}

/**
   @brief Evaluates the status : charging, discharging or idle phase of the storage medium
   @param[in] Time evaluated from the starting of the cycle (incremented by time_step)
   @param[out] 1 - Charging phase, 2 - Discharging phase, 3 - Idle phase
**/

int status_storage(double c)
{		
	double time_discharge_start, time_discharge_end;	
	time_discharge_start = time_charge + time_idle;
	time_discharge_end = time_discharge_start + time_discharge;
	
	if(c <= time_charge)
		return 1; //status - charging (1)
	else if(c > time_discharge_start && c < time_discharge_end)
		return 2; //status - discharging (2)
	else
		return 3; //status - idle (3)
	
}

/**
   @brief Solves the discretized equations of the fluid - both MMS and without
   @param[in] Status of the phase (Charging/Discharging/Idle)
   @param[in] Wavenumber specifically chosen for the fluid - k_waveno1
**/

void liq_solve(int a, double k_waveno)
{
  double s_0,s_n;

  if(ovs_liquid == 1) //MMS - only charging phase
    {
      //computing the source terms at the boundaries
      s_0 = - (u_f*k_waveno*sin(k_waveno*0.5*delta_x)) + (k_waveno*k_waveno*diff_liq*cos(k_waveno*0.5*delta_x));
      s_n = - (u_f*k_waveno*sin(k_waveno*((double)(grid_no-1)*delta_x + 0.5*delta_x))) + 
	(k_waveno*k_waveno*diff_liq*cos(k_waveno*((double)(grid_no-1)*delta_x + 0.5*delta_x)));
     
      t_fluid_new[0] = t_fluid[0] - (u_f*time_step*(t_fluid[0] - 1))/delta_x + 
	(diff_liq*time_step*(t_fluid[1] - t_fluid[0])/(delta_x*delta_x)) + s_0*time_step;
     
      for(int i=1; i < grid_no-1; i++)
	  t_fluid_new[i] = t_fluid[i] - (u_f*time_step*(t_fluid[i] - t_fluid[i-1]))/delta_x  + 
	    (diff_liq*time_step*(t_fluid[i+1] - 2*t_fluid[i] + t_fluid[i-1]))/(delta_x*delta_x) - 
	    (u_f*time_step*k_waveno*sin(k_waveno*0.5*(double)(2*i + 1)*delta_x)) + 
	    (k_waveno*k_waveno*diff_liq*time_step*cos(k_waveno*0.5*(double)(2*i + 1)*delta_x));
		
      t_fluid_new[grid_no - 1] = t_fluid[grid_no - 1] - (u_f*time_step*(t_fluid[grid_no-1] - t_fluid[grid_no-2]))/delta_x
		  + (diff_liq*time_step*(t_fluid[grid_no - 2] - t_fluid[grid_no - 1])/(delta_x*delta_x)) + s_n*time_step;
    }

  else
    {
      if(a == 1) //Charging phase - normal simulation
	{      
	  t_fluid_new[0] = t_fluid[0] - (u_f*time_step*(t_fluid[0] - t_fluid_bc_left))/delta_x + 
	    (diff_liq*time_step*(t_fluid[1] - t_fluid[0])/(delta_x*delta_x));

	  for(int i=1; i < grid_no-1; i++)
	    t_fluid_new[i] = t_fluid[i] - (u_f*time_step*(t_fluid[i] - t_fluid[i-1]))/delta_x  + 
	      (diff_liq*time_step*(t_fluid[i+1] - 2*t_fluid[i] + t_fluid[i-1]))/(delta_x*delta_x) ;
		
	  t_fluid_new[grid_no - 1] = t_fluid[grid_no - 1] - (u_f*time_step*(t_fluid[grid_no-1] - t_fluid[grid_no-2]))/delta_x
	    + (diff_liq*time_step*(t_fluid[grid_no - 2] - t_fluid[grid_no - 1])/(delta_x*delta_x));
	}
  
      else if(a == 2) //Discharging phase - normal simulation
	{
	  t_fluid_new[0] = t_fluid[0] + (u_f*time_step*(t_fluid[1] - t_fluid[0]))/delta_x
	    + (diff_liq*time_step*(t_fluid[1] - t_fluid[0])/(delta_x*delta_x));
	  
	  for(int i=1; i < grid_no-1; i++)
	    t_fluid_new[i] = t_fluid[i] + (u_f*time_step*(t_fluid[i+1] - t_fluid[i]))/delta_x 
	      + (diff_liq*time_step*(t_fluid[i+1] - 2*t_fluid[i] + t_fluid[i-1]))/(delta_x*delta_x);
	        
	  t_fluid_new[grid_no - 1] = t_fluid[grid_no - 1]
	    + (u_f*time_step*(t_fluid_bc_right - t_fluid[grid_no-1]))/delta_x  
	    + (diff_liq*time_step*(t_fluid[grid_no - 2] - t_fluid[grid_no - 1]))/(delta_x*delta_x);
	}
			
      else  //Idle phase - normal simulation
	{

	  t_fluid_new[0] = t_fluid[0]  + (diff_liq*time_step*(t_fluid[1] - t_fluid[0])/(delta_x*delta_x));
     
          for(int i=1; i < grid_no-1; i++)
	    t_fluid_new[i] = t_fluid[i]  + 
	      (diff_liq*time_step*(t_fluid[i+1] - 2*t_fluid[i] + t_fluid[i-1]))/(delta_x*delta_x);
		
	  t_fluid_new[grid_no - 1] = t_fluid[grid_no - 1]
		  + (diff_liq*time_step*(t_fluid[grid_no - 2] - t_fluid[grid_no - 1])/(delta_x*delta_x)) ;
	} 
    }
}

/**
   @brief Solves the coupled solid and fluid equations - point implicit method
   @param[in] Status of the phase (Charging/Discharging/Idle)
   @param[in] Wavenumber specifically chosen for the fluid - k_waveno1
   @param[in] Wavenumber specifically chosen for the solid - k_waveno2
**/

void coupled_solve_implicit(int a, double k_waveno1, double k_waveno2)
{

  double diff;
  liq_solve(a, k_waveno1); //solves the fluid eqns and updates the t_fluid_new array
  solid_solve(a, k_waveno2); //solves the solid eqns and updates the t_solid_new array

  if(ovs_liquid == 1 && ovs_solid == 1) //if MMS for both - Part 4
    {
      for(int i = 0; i < grid_no; ++i)
	{
	  //correcting the source terms by incorporating them here, as thse have not been accounted in liq_solve() and solid_solve()
	  diff = cos(k_waveno1*0.5*(double)(2*i + 1)*delta_x) - cos(k_waveno2*0.5*(double)(2*i + 1)*delta_x);
	  t_fluid_new[i] = t_fluid_new[i] + (diff*h_vf*time_step);
	  t_solid_new[i] = t_solid_new[i] - (diff*h_vs*time_step); //negative sign as negative of diff has to be taken : cosk2 - cosk1
	}
    }

  //Solves the 2 equations at each point using Gaussian Elimination
      for(int i = 0; i < grid_no; ++i)  
	  gauss_elimination(&(t_fluid_new[i]),&(t_solid_new[i]));
	
 }
	
/**
   @brief Solves the discretized equations of the solid - both MMS and without
   @param[in] Status of the phase (Charging/Discharging/Idle)
   @param[in] Wavenumber specifically chosen for the solid - k_waveno2
**/

void solid_solve(int a, double k_waveno)
{	
  double s_0,s_n;
  
  if(ovs_solid == 1) //MMS - charging phase
    {
      //computing the source terms at the boundaries
      s_0 =  k_waveno*k_waveno*diff_solid*cos(k_waveno*0.5*delta_x);
      s_n =  k_waveno*k_waveno*diff_solid*cos(k_waveno*((double)(grid_no-1)*delta_x + 0.5*delta_x));
     
     t_solid_new[0] = t_solid[0] + (diff_solid*time_step*(t_solid[1] - t_solid[0])/(delta_x*delta_x)) 
       + s_0*time_step;
     
     for(int i=1; i < grid_no-1; i++)
       t_solid_new[i] = t_solid[i] 
	 + (diff_solid*time_step*(t_solid[i+1] - 2*t_solid[i] + t_solid[i-1]))/(delta_x*delta_x)  
	 + (k_waveno*k_waveno*diff_solid*time_step*cos(k_waveno*0.5*(double)(2*i + 1)*delta_x));
     
     t_solid_new[grid_no - 1] = t_solid[grid_no - 1]  
       + (diff_solid*time_step*(t_solid[grid_no - 2] - t_solid[grid_no - 1]))/(delta_x*delta_x) 
       + s_n*time_step;
    }

  else  // normal simulation 
    {      
     t_solid_new[0] = t_solid[0] + (diff_solid*time_step*(t_solid[1] - t_solid[0])/(delta_x*delta_x));
     
     for(int i=1; i < grid_no-1; i++)
       t_solid_new[i] = t_solid[i] + (diff_solid*time_step*(t_solid[i+1] - 2*t_solid[i] 
          + t_solid[i-1]))/(delta_x*delta_x) ;
     
     t_solid_new[grid_no - 1] = t_solid[grid_no - 1]  
       + (diff_solid*time_step*(t_solid[grid_no - 2] - t_solid[grid_no - 1]))/(delta_x*delta_x) ;
    }

}

/**
   @brief Computes the error for OVS - solid
   @param[in] Wavenumber specifically of the solid - k_waveno2
**/

void error_compute_ovs_solid(double k_waveno)
{
  double error_l1 = 0.0,error_l2 = 0.0, error_linf = 0, term;
  int pos = 0;
  FILE* fileID, *file_temp;
  fileID = fopen("error_ovs_solid.dat","a"); //file with the l1, l2, l-inf error
  file_temp = fopen("err_solid_var.dat","w"); //file with the temp diff data

  for(int j = 0; j < grid_no; j++)
	{
	        term = t_solid[j] - cos(k_waveno*delta_x*0.5*(double)(2*j + 1));
	        fprintf(file_temp,"%d\t %f\n",j, fabs(term));
		error_l1 += fabs(term);
		error_l2 += term*term;
		if(fabs(term) > error_linf)
		{
			error_linf = fabs(term);		
				pos = j;
		}
	}

        error_l1 = error_l1/(double)grid_no;
        error_l2 = sqrt(error_l2/(double)grid_no);
	fprintf(fileID,"%lf\t %7.9e\t %7.9e\t %7.9e\t %d\n",delta_x,error_l1,error_l2,error_linf,pos);
	fclose(fileID);
	fclose(file_temp);
}

/**
   @brief Computes the error for OVS - fluid
   @param[in] Wavenumber specifically of the fluid - k_waveno1
**/ 

void error_compute_ovs_liquid(double k_waveno)
{
        double error_l1 = 0.0,error_l2 = 0.0, error_linf = 0, term;
	int pos = 0; //required for noting the location of maximum error 
	FILE* fileID, *file_temp;

	fileID = fopen("error_ovs_fluid.dat","a"); //file with the l1, l2, l-inf error
	file_temp = fopen("err_fluid_var.dat","w"); //file with the temp diff data 
	
	//Computing the errors for fluid 
	for(int j = 0; j < grid_no; j++)
	{
	        term = t_fluid[j] - cos(k_waveno*delta_x*0.5*(double)(2*j + 1));
	        fprintf(file_temp,"%d\t %f\n",j, fabs(term));
		error_l1 += fabs(term);
		error_l2 += term*term;
		if(fabs(term) > error_linf)
		{
			error_linf = fabs(term);		
			pos = j;
		}
	}
       
       error_l1 = error_l1/(double)grid_no;
       error_l2 = sqrt(error_l2/(double)grid_no);
       fprintf(fileID,"%lf\t %7.9e\t %7.9e\t %7.9e\t %d\n",delta_x,error_l1,error_l2,error_linf, pos); 
	
       fclose(fileID);
       fclose(file_temp);
	
}

/**
   @brief Computes the error to note the convergence - required for OVS 
   @param[out] Error(Average) in double precision
**/

double conv_error()
{
                double error_mean_1 = 0.0, error_mean_2 = 0.0;
		for(int j = 0; j < grid_no; j++)
		  {
		    	error_mean_1 += fabs(t_fluid_iter[j] - t_fluid[j]);
			error_mean_2 += fabs(t_solid_iter[j] - t_solid[j]);
		  }
		error_mean_1 = error_mean_1/(double)grid_no; //averaged over all the grid cells
		error_mean_2 = error_mean_2/(double)grid_no;

		if(error_mean_1 < error_mean_2)
		  error_mean_1 = error_mean_2; //chooses the maximum error 
		
		return error_mean_1; //largest error is returned
}
  
/**
   @brief Computes the maximum amount of thermal energy stored 
   @param[out] Thermal Energy maximum amount in double precision
**/

double q_max()
{
  
  double vol, max_energy;
  vol = M_PI*diameter*diameter*height/4.0;
  max_energy = (epsilon*rho_f*C_f + (1.0 - epsilon)*rho_s*C_s)*vol*
    (t_fluid_bc_left - t_fluid_bc_right);
  
  return max_energy;
}
  

/**
   @brief Computes the amount of thermal energy stored at a certain time instant
   @param[out] Thermal Energy amount in double precision
**/
	       
double q_time()
{
  double area, solid_sum = 0., liq_sum = 0., q;
  area = M_PI*diameter*diameter/4.0;
  
  for(int i = 0; i < grid_no; ++i)
    {
      liq_sum += t_fluid[i] - t_fluid_bc_right;
      solid_sum += t_solid[i] - t_fluid_bc_right;
    }
  
  q = (epsilon*rho_f*C_f*liq_sum*delta_x + (1 - epsilon)*rho_s*C_s*solid_sum*delta_x)*area;
  
  return q;
}

/**
   @brief Computes the exergy flux
   @param[in] Position at which this exergy flux has to be caclulated
   @param[out] Exergy Flux in double precision
**/

double exergy_flux(int pos)
{
  double exergy;
  exergy = mass_flow*C_f*(t_fluid[pos] - 288.15 - 288.15*log(t_fluid[pos]/288.15))*time_step;

  return exergy;
}
  
  
