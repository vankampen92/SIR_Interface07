/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                        MEASLES PARAMETERS                           */
/*                          (PNAS, 2004, 101:16915-16916)                    */
/*                             CONSTANT COMMUNITY SIZE                       */
/*                                                                           */
/*                             David Alonso, 2005 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../SIR.h"
#include "../../SIR_Analytic_General.h"
int No_of_Points;
/* Global Shared parameters main Program <---> ArgumentControl() */
int TIMES;     /* Number of times to be analyzed */
float STEP_SIZE; 
int Realizations;
int Simlength;
int I_0, M_0;      /* Inital number of infective Individuals */
double a,d,g,mu,b,Imm;
int POPULATION;
double Per;  /* Period of the seasonal Forcing: usually 365 days */
double b_0, b_1;
float l_0, l_1; /* Range in Betas (transmission rates) */
float d_0, d_1; /* Range in Deltas (turnover rates) */
float r_0, r_1; /* Range in Immigration rates */
int TIME_DEPENDENT_PARAM;
float time_Factor; /* Time (dimless units) = time (days) * time_Factor; */
/* End Global Shared parameters main Program <---> ArgumentControl() */

int no_Main_Max;
void Initial_Settings_Measles();

void Initial_Settings_Measles()
{
  TIME_DEPENDENT_PARAM = 0; no_Main_Max = 10;
  Realizations = 2000; No_of_Points =  2000;
  Simlength = 8000000;
  POPULATION = 500000;
  TIMES = 500; STEP_SIZE = 14.;
  time_Factor = 1./13.;
  b = 15.27; a = 0.; d = 7.15e-4; g = 1.; mu = 0.; Imm = 1.3e-4;
  b_0 = b; b_1 = 0.25; Per = 365.*time_Factor;
  I_0 = 1380; M_0 = 447970;
  /* Dimensionless Parameters */
  l_0 = 10.;                   l_1 = 20.;               /* Range in Betas (transmission rates) */
  d_0 = 6.0e-4;                d_1= 2.0e-3;             /* Range in Deltas (turnover rates) */
  r_0 = 1.e-5/time_Factor;     r_1 = 2.e-5/time_Factor; /* Range in Immigration rates */
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
  float *x,*y;
  float Mean, Var;

  //Initial_Settings_DushoffEtal();  
  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  Initial_Settings_Measles();
  InitSeed();
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */
  
  /* Command line arguments */
  printf(" Parameters should be introduced in dimensionless units\n");
  if(argc>1) ArgumentControl(argc,argv);

  range_L[0] = l_0; range_L[1] = l_1;
  range_D[0] = d_0; range_D[1] = d_1;
  range_R[0] = r_0; range_R[1] = r_1;
  
  x = vector(1,No_of_Points); y = vector(1,No_of_Points);
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
    sample_Parameters(&Par, range_L, range_D, range_R);
    Fixed_Points_General(&Par, &f_0, &f_1); Stability_General(&Par);
    /* Initial Equilibrium State */
    S_0 = f_0 * POPULATION; I_0 = f_1 * POPULATION; 
    printf("Initial Equilibrium State: S = %d\tI = %d\n", S_0, I_0);
  
    M_0 = POPULATION - S_0 - I_0;
    printf("Transmission, b_0 = %g\tTurnover rate (death rate) = %g\tImmigration level= %g\n",
	   b_0, Par.Delta, Par.Imm);
    period[i] = 1./time_Factor /Resonance_Frequency_Peak(&Par, 1) / 356.;
    
    /* Seasonal Focing */
    TIME_DEPENDENT_PARAM = 3;
    printf("School-term seasonal forcing: S = %d\n", TIME_DEPENDENT_PARAM);
    SIR_StocRea(Village, &Par, Simlength, modul, 0., S_I, Time, Yp, I_0);

    for(j=0; j<No_of_Points; j++) y[j+1] = (float)Yp[j][1];
  
    Mean = mean_amplitude(y,No_of_Points);
    Var  = mean_squared_amplitude(y,No_of_Points) - Mean*Mean;

    //ratio[i] = sqrt(Var)/Mean;
    ratio[i] = sqrt(Var);
    
    printf("Final State: Time = %g\tNumber of Infectious Individuals = %d\n", 
	   Time[No_of_Points-1], Yp[No_of_Points-1][1]);
    printf("Endogenous Period: %f\t Sqrt(Var) (Seasonal Forcing): %f\n\n", 
	   period[i], ratio[i]);
    
    /* Non Seasonal Forcing for the same parameter choice */
    TIME_DEPENDENT_PARAM = 0; 
    printf("Non-seasonal forcing: S = %d\n", TIME_DEPENDENT_PARAM);
    SIR_StocRea(Village, &Par, Simlength, modul, 0., S_I, Time, Yp, I_0);
    
    for(j=0; j<No_of_Points; j++) y[j+1] = (float)Yp[j][1];
  
    Mean = mean_amplitude(y,No_of_Points);
    Var  = mean_squared_amplitude(y,No_of_Points) - Mean*Mean;

    //RATIO[i] = sqrt(Var)/Mean;
    RATIO[i] = sqrt(Var);
    
    printf("Final State: Time = %g\tNumber of Infectious Individuals = %d\n", 
	   Time[No_of_Points-1], Yp[No_of_Points-1][1]);
    printf("Endogenous Period: %f\t Sqrt(Var) (Non_Forcing): %f\n", 
	   period[i], RATIO[i]); 

    fprintf(fp,"%g\t%g\t%g\t%f\n", b_0, Par.Delta, Par.Imm, ratio[i]);
    fprintf(fp0,"%g\t%g\t%g\t%f\n", b_0, Par.Delta, Par.Imm, RATIO[i]);
    printf("... ... ... ... ... ... ... ...\n");
  }
  fclose(fp); fclose(fp0);
  /* End of STOCHASTIC REALIZATIONS for one particular time * * * * * * * */ 
  Saving_to_File_float_float("rat_", period, ratio, Realizations,1,0);
  Saving_to_File_float_float("rat_", period, RATIO, Realizations,0,0);
  Saving_to_File_float_float("r.TO.R", RATIO, ratio, Realizations,0,0); 
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  free_imatrix(Yp, 0,No_of_Points, 0, 2);
  free_dvector(Time, 0, No_of_Points);
  free_dvector(Time_Max, 0,No_of_Points); free_dvector(Time_Min, 0,No_of_Points);
  free_dvector(Maxim, 0,No_of_Points); free_dvector(Minim, 0,No_of_Points);
  free_vector(ratio, 0,Realizations); free_vector(period, 0,Realizations);
  free_vector(RATIO, 0,Realizations); 
  
  free_vector(x, 1,No_of_Points); free_vector(y, 1,No_of_Points);
  printf("\n\n\n\nNon-seasonal forcing: S = %d\n", TIME_DEPENDENT_PARAM);
  modelReport("report.txt");

  printf("\nEnd of progam\n");
  return (0);
}







