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
extern int TIME_DEPENDENT_PARAM;

float det_season_amplif(int MSD, ParamSet *Pa, float T_sup, int no_of_Maxima, int SI)
{
  int N,i;              /* Number of equations      */
  int No_of_Points;  /* Number of discret points */ 
  float p_1, p_2, x_i, x_s;
  double Fi_0, Psi_0;
  float *P, *x, *y;
  double *Time,time;
  double **X; int **Yp;
  double *Time_Max, *Time_Min, *Maxim, *Minim;
  int no_Max, I_Level_0;
  float ratio, Var, Mean;

  N = 2; /* Number of Equations */
  POPULATION = 100000;     
  No_of_Points = 1000;  /* Number of integration points */

  b = Pa->Beta; a = Pa->Alpha; d = Pa->Delta; g = Pa->Gamma; mu = Pa->Mu; 
  Imm = Pa->Imm;  b_0 = Pa->Beta; b_1 = 0.25;
  x_i = 0.; x_s = T_sup;
  
  Time_Max = dvector(0,No_of_Points); Time_Min = dvector(0,No_of_Points);
  Maxim = dvector(0,No_of_Points); Minim = dvector(0, No_of_Points);
  Time = dvector(0,No_of_Points); Yp = imatrix(0,No_of_Points, 0, N-1);
  X = dmatrix(0,No_of_Points, 0,N-1);
  x = vector(1,No_of_Points); y = vector(1,No_of_Points); 
  P = vector(0,N-1);
  //Fixed_Points_General(Pa, &Fi, &Psi);
  Fi_0 = (float)Fi(Pa); Psi_0 = (float)Psi(Pa);
  P[0] = Fi_0; P[1] = Psi_0; 
   
  integration_01(P, N, x_i, x_s, No_of_Points, deriva, 
		 Time, X);
  if(MSD == 0){
    for(i=0; i<No_of_Points; i++){
      Yp[i][0] = (int)(POPULATION*X[i][0]); 
      Yp[i][1] = (int)(POPULATION*X[i][1]); 
    } 
    
    I_Level_0 = (int)(POPULATION*Psi_0);
    set_to_value_double(Maxim,No_of_Points,0.);
    set_to_value_double(Minim,No_of_Points,0.);
    peak_to_Trough_Maxims_Minims(Time, Yp, No_of_Points, I_Level_0,
				 Time_Max, Time_Min, Maxim, Minim, &no_Max);
    if(no_Max > no_of_Maxima)
      the_Largest_Maxima(Maxim, no_Max, no_of_Maxima);
    else
      no_of_Maxima = no_Max;
    /* Average amplification calculated with a fixed $no_Max$ peak points */
    ratio = Relative_Amplification_Ratio(Maxim, no_of_Maxima, I_Level_0); 
  }
  
  else{
    for(i=0; i<No_of_Points; i++){
      x[i+1] = X[i][0]; 
      y[i+1] = X[i][1]; 
    } 
    
    if(SI == 0){
      Mean = mean_amplitude(x,No_of_Points);
      Var = mean_squared_amplitude(x,No_of_Points) - Mean*Mean; 
    }
    else if(SI == 1){
      Mean = mean_amplitude(y,No_of_Points);
      Var = mean_squared_amplitude(y,No_of_Points) - Mean*Mean;
    }
    else{
      printf(" The system is two dimensional: either SI = 0 or SI = 1, but SI = %d!!!\n", 
	     SI);
      exit(0);
    }
    ratio = sqrt(Var)/Mean;
  }

  free_vector(x,1,No_of_Points); free_vector(y, 1, No_of_Points);
  free_imatrix(Yp, 0,No_of_Points, 0, N-1);
  free_dmatrix(X, 0,No_of_Points, 0,N-1);
  free_dvector(Time, 0, No_of_Points);
  free_dvector(Time_Max, 0,No_of_Points); free_dvector(Time_Min, 0,No_of_Points);
  free_dvector(Maxim, 0,No_of_Points); free_dvector(Minim, 0,No_of_Points);
  free_vector(P,0,N-1);
  return (ratio);
}
