/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*                      Computing Deterministic Amplification                */
/*                             David Alonso, 2005 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../SIR.h"
#include "../../SIR_Analytic_General.h"

extern double a,d,g,mu,b,Imm;
extern double b_0,b_1;
extern int POPULATION;
int I_0, M_0;
float t_0;

float sto_season_amplif(int MSD, ParamSet *Pa, float T_sup, int no_of_Maxima, int SI)
{
  int N, i, j, new;              /* Number of equations      */
  int No_of_Points;  /* Number of discret points */ 
  float x_i, x_s;
  float *x, *y;
  double Fi_0, Psi_0;
  double STEP_SIZE;
  double *Time;
  int **Yp;
  double *Time_values, *Time_Max, *Time_Min, *Maxim, *Minim;
  int no_Max, I_Level_0;
  float ratio, Var, Mean;
  int Yp_End[1][3], Strain[3]; /*0: S; 1: I; 2: R */
  double Time_End[1];
  Community Village[No_of_Villages];

  N = 3; /* Number of Different Species */    
  No_of_Points = 500;  /* Number of discretization time points */

  b = Pa->Beta; a = Pa->Alpha; d = Pa->Delta; g = Pa->Gamma; mu = Pa->Mu; 
  Imm = Pa->Imm;  b_0 = Pa->Beta; b_1 = 0.25;
  t_0 = x_i = 0.; x_s = T_sup;
  
  Time_values = dvector(0,No_of_Points-1);
  Time_Max = dvector(0,No_of_Points); Time_Min = dvector(0,No_of_Points);
  Maxim = dvector(0,No_of_Points); Minim = dvector(0, No_of_Points);
  Time = dvector(0,No_of_Points); Yp = imatrix(0,No_of_Points, 0, N-1);
  x = vector(1,No_of_Points); y = vector(1,No_of_Points); 

  /* Fixed Poitns */
  Fi_0 = (float)Fi(Pa); /* S */ Psi_0 = (float)Psi(Pa); /* I */
  I_0 = POPULATION*Psi_0; M_0 = POPULATION*(1.- Fi_0 - Psi_0);

  STEP_SIZE = ((double)T_sup - (double)x_i)/(double)No_of_Points;
  /*setting times*/
  for(i=0; i<No_of_Points; i++)
    Time_values[i] = x_i + (double)(i+1) * STEP_SIZE;

  //printf("Initial population values: S = %d\t I = %d\n", POPULATION -I_0-M_0, I_0);
  //printf("Sampling Interval (days): %g\n", STEP_SIZE);
  //printf("Entering Generation of Stochastic Realizations...\n");
  //Press_Key();
    
  Time_End[0] = x_i;
  for (j=0; j<No_of_Points; j++){

    SIR(Village, Pa, 1, 1, Time_values[j], Strain, Time_End, Yp_End, &new);
    
    Time[j] = Time_End[0];  Yp[j][0] = Strain[0]; Yp[j][1] = Strain[1];
    
    /* Time_End[0] is always the last time which is the closest to Time_values[j] */
    /* ... ... going further for the next time */
  }

  if(MSD == 0){
    I_Level_0 = (int)(POPULATION*Psi_0);
   set_to_value_double(Maxim,No_of_Points,0.);
   set_to_value_double(Minim,No_of_Points,0.); 
   peak_to_Trough_Maxims_Minims(Time, Yp, No_of_Points, I_Level_0,
				 Time_Max, Time_Min, Maxim, Minim, &no_Max);
    the_Largest_Maxima(Maxim, no_Max, no_of_Maxima);
    /* Average amplification calculated with a fixed $no_of_Maxima$ peak points */
    ratio = Relative_Amplification_Ratio(Maxim, no_of_Maxima, I_Level_0); 
  }
  else{
    for(i=0; i<No_of_Points; i++){
      x[i+1] = (float)Yp[i][0]; 
      y[i+1] = (float)Yp[i][1]; 
    } 
    
    if(SI == 0){
      Mean = mean_amplitude(x,No_of_Points);
      Var = mean_squared_amplitude(x,No_of_Points) - 2.*(float)POPULATION*Fi_0*Mean + (float)POPULATION*Fi_0*(float)POPULATION*Fi_0; 
    }
    else if(SI == 1){
      Mean = mean_amplitude(y,No_of_Points);
      Var = mean_squared_amplitude(y,No_of_Points) - 2.*(float)POPULATION*Psi_0*Mean+ (float)POPULATION*Psi_0*(float)POPULATION*Psi_0;
    }
    else{
      printf(" The system is two dimensional: either SI = 0 or SI = 1, but SI = %d!!!\n", 
	     SI);
      exit(0);
    }
    /* Average amplification calculated as the coeficient of variantion of the time series */
    ratio = Var/(double)POPULATION;
  }

  free_vector(x,1,No_of_Points); free_vector(y, 1, No_of_Points);
  free_imatrix(Yp, 0,No_of_Points, 0, N-1);
  free_dvector(Time, 0, No_of_Points);
  free_dvector(Time_Max, 0,No_of_Points); free_dvector(Time_Min, 0,No_of_Points);
  free_dvector(Maxim, 0,No_of_Points); free_dvector(Minim, 0,No_of_Points);
  free_dvector(Time_values, 0,No_of_Points-1);
  return (ratio);
}
/* Writing amplification...*/
//for(i=0; i<no_Max; i++) printf("Max[%d] = %g, ", i, Maxim[i]);
//printf("\n");
//Press_Key();
