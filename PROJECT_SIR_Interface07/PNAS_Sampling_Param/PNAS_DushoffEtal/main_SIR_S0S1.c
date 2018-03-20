/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (PNAS, 2004, 101:16915-16916)                    */
/*                             CONSTANT COMMUNITY SIZE                       */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../SIR.h"
int No_of_Points;
/* Global Shared parameters main Program <---> ArgumentControl() */
int TIMES;     /* Number of times to be analyzed */
float STEP_SIZE; 
float EPSILON;
int Realizations;
int Simlength;
int I_0, M_0;      /* Inital number of infective Individuals */
double a,d,g,mu,b,Imm;
int POPULATION;
double Per;  /* Period of the seasonal Forcing: usually 365 days */
double b_0, b_1;
float l_0, l_1;
float d_0, d_1;
float r_0, r_1;
float sensibility;
int SENSIBILITY;
int TIME_DEPENDENT_PARAM;
/* End Global Shared parameters main Program <---> ArgumentControl() */

float time_Factor;
int no_Main_Max;
void Initial_Settings_DushoffEtal()
{
  TIME_DEPENDENT_PARAM = 0; SENSIBILITY = 0; sensibility = 0.04; no_Main_Max = 10;
  Realizations = 2000; No_of_Points =  2000; time_Factor = 1.;
  Simlength = 8000000;
  POPULATION = 500000;
  TIMES = 500; STEP_SIZE = 14.; EPSILON = 0.001;
  
  b = 400.; a = 0.; d = 0.; g = 1./0.025; mu = 1/8.; Imm = 0.;
  b_0 = 400.; b_1 = 0.04; Per = 1.;
  I_0 = 1380; M_0 = 447970;
  l_0 = 4.;       l_1 = 8.;
  d_0 = 6.0/365.; d_1= 10.0/365.;
  r_0 = 4.;       r_1 = 16.;
}

int main(int argc, char **argv)
{
  int N, S_0; /* Total Metapopulation Capacity */
  int i,j,k,modul, count, no_Max; 
  int S_I[3];
  int **Yp;
  double *Time;
  double f_0, f_1;
  float DELTA;
  float *ratio, *period, *RATIO;
  float range_D[2], range_L[2], range_R[2];
  /* Temporal evolution of the number of strains in each stochastic realization */
  Community Village[No_of_Villages];
  ParamSet Par;
  FILE *fp, *fp0;
  double *Time_Max, *Time_Min, *Maxim, *Minim;

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  Initial_Settings_DushoffEtal();
  InitSeed();
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

  range_L[0] = l_0; range_L[1] = l_1;
  range_D[0] = d_0; range_D[1] = d_1;
  range_R[0] = r_0; range_R[1] = r_1;
  Yp = imatrix(0,No_of_Points, 0, 2);
  Time = dvector(0, No_of_Points);
  Time_Max = dvector(0,No_of_Points); Time_Min = dvector(0,No_of_Points);
  Maxim = dvector(0,No_of_Points); Minim = dvector(0, No_of_Points);
  ratio = vector(0, Realizations); period = vector(0, Realizations);
  RATIO = vector(0, Realizations);  
  N = POPULATION; /* Village Size */ 
  settingParameterStruct(&Par);
  modelReport("report.txt");
  
  if(Simlength < No_of_Points){
    printf("Simulation length must be at least equal to the number of points to be saved\n");
    exit(0);
  }
  if(Simlength%No_of_Points == 0)
    modul = Simlength/No_of_Points;
  else
    modul = Simlength/No_of_Points + 1;
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* Intiating STOCHASTIC REALIZATIONS */
  fp = fopen("parameter1.dat", "w"); fp0 = fopen("parameter0.dat", "w");
  for (i=0;i<Realizations;i++){
    printf("REALIZATION: %dth!!!!\n", i); 
    sample_Parameters(&Par, range_R, range_L, range_D);
    Fixed_Points_General(&Par, &f_0, &f_1); Stability_General(&Par);
    /* Initial Equilibrium State */
    S_0 = f_0 * POPULATION; I_0 = f_1 * POPULATION; 

    if(SENSIBILITY == 1){
      printf("Deterministic Equilibrium State: S = %d\tI = %d\n", S_0, I_0);
      if(drand48() > 0.5) DELTA = +sensibility * (float)I_0;
      else                DELTA = -sensibility * (float)I_0;
      I_0 += DELTA; S_0 -= DELTA; 
      printf("Initial Fuzzy Equilibrium State: S = %d\tI = %d\n", S_0, I_0);
    }
    else printf("Initial Equilibrium State: S = %d\tI = %d\n", S_0, I_0);
  
    M_0 = POPULATION - S_0 - I_0;
    printf("Transmission, b_0 = %g\tInfectious Period, D = %g\tImmunity Period = %g\n",
	   b_0, 365./Par.Gamma, 1./Par.Mu);
    period[i] = approximate_Endogenous_Period(&Par);

    /* Seasonal Focing */
    TIME_DEPENDENT_PARAM = 1;
    printf("School term seasonal forcing: S=%d\n", TIME_DEPENDENT_PARAM);
    SIR_StocRea(Village, &Par, Simlength, modul, 0., S_I, Time, Yp, I_0);
    if(Yp[No_of_Points-1][1] > 0){
      set_to_value_double(Maxim,No_of_Points,0.);
      set_to_value_double(Minim,No_of_Points,0.);
      peak_to_Trough_Maxims_Minims(Time, Yp, No_of_Points, I_0,
				   Time_Max, Time_Min, 
				   Maxim, Minim, &no_Max);

      no_Main_Max = 10; 
      ordering_Vector_2_double(Maxim, Time_Max, no_Max, -1);
      ordering_Vector_2_double(Time_Max, Maxim, no_Main_Max, 1);
      if(SENSIBILITY == 0)
	ratio[i] = Relative_Amplification_Ratio_Forcing(Time_Max, Maxim, no_Main_Max, 
							POPULATION, f_0,f_1, &Par, 
							No_of_Points);
    }
    else ratio[i] = 0.;
    printf("Final State: Time = %g\tNumber of Infectious Individuals = %d\n", 
	   Time[No_of_Points-1], Yp[No_of_Points-1][1]);
    printf("Endogenous Period: %f\tPeak to Trough Ratio (Seasonal Forcing): %f\n", 
	   period[i], ratio[i]); 
    /* Non Seasonal Forcing */
    TIME_DEPENDENT_PARAM = 0;  
    printf("Non-seasonal forcing: S=%d\n", TIME_DEPENDENT_PARAM);
    SIR_StocRea(Village, &Par, Simlength, modul, 0., S_I, Time, Yp, I_0);
    if(Yp[No_of_Points-1][1] > 0){
      set_to_value_double(Maxim,No_of_Points,0.);
      set_to_value_double(Minim,No_of_Points,0.);
      peak_to_Trough_Maxims_Minims(Time, Yp, No_of_Points, I_0,
				   Time_Max, Time_Min, 
				   Maxim, Minim, &no_Max);
      
      no_Main_Max = 10; ordering_Vector_double(Maxim, no_Max, -1);
      RATIO[i] = Relative_Amplification_Ratio(Maxim, no_Main_Max, I_0);
    }
    else RATIO[i] = 0.;
    printf("Final State: Time = %g\tNumber of Infectious Individuals = %d\n", 
	   Time[No_of_Points-1], Yp[No_of_Points-1][1]);
    printf("Endogenous Period: %f\tPeak to Trough Ratio (Non_Forcing): %f\n", 
	   period[i], RATIO[i]); 

    fprintf(fp,"%g\t%g\t%g\t%f\n", b_0, Par.Gamma, Par.Mu, ratio[i]);
    fprintf(fp0,"%g\t%g\t%g\t%f\n", b_0, Par.Gamma, Par.Mu, RATIO[i]);
    printf("... ... ... ... ... ... ... ...\n");
  }
  fclose(fp); fclose(fp0);
  /* End of STOCHASTIC REALIZATIONS for one particular time * * * * * * * */ 
  Saving_to_File_float_float("rat_", period, ratio, Realizations,1,0);
  Saving_to_File_float_float("rat_", period, RATIO, Realizations,0,0);
  Saving_to_File_float_float("rTOR", RATIO, ratio, Realizations,0,0); 
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  free_imatrix(Yp, 0,No_of_Points, 0, 2);
  free_dvector(Time, 0, No_of_Points);
  free_dvector(Time_Max, 0,No_of_Points); free_dvector(Time_Min, 0,No_of_Points);
  free_dvector(Maxim, 0,No_of_Points); free_dvector(Minim, 0,No_of_Points);
  free_vector(ratio, 0,Realizations); free_vector(period, 0,Realizations);
  free_vector(RATIO, 0,Realizations); 
  modelReport("report.txt");

  printf("\nEnd of progam\n");
  return (0);
}







