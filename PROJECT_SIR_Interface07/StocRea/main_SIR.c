/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../SIR.h"
int No_of_Points;
/* Global Shared parameters main Program <---> ArgumentControl() */
int Realizations;
int Simlength;
int I_0, M_0;      /* Inital number of infective Individuals */
double a,d,g,mu,b,Imm;
int POPULATION;
double Per;  /* Period of the seasonal Forcing: usually 365 days */
double b_0, b_1;
int TIME_DEPENDENT_PARAM;
float time_Factor;

int main(int argc, char **argv)
{
  int N; /* Total Metapopulation Capacity */
  int i,j,k,modul, count, spNumber; 
  int S_I[3];
  int **Yp;
  double *Time;
  double f_0, f_1;
  /* Temporal evolution of the number of strains in each stochastic realization */
  Community Village[No_of_Villages];
  ParamSet Par;
  double *Time_Max, *Time_Min, *Maxim, *Minim;

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  TIME_DEPENDENT_PARAM = 0;
  Realizations = 1; No_of_Points =  2000;
  Simlength = 500000;
  POPULATION = 100000;
  time_Factor = 1./13.;	
  InitSeed();
  b = 18.; a = 0.; d = 1.64e-4; g = 1.; mu = 0.; Imm = 1.98e-4;
  b_0 = b; b_1 = 0.25; Per = 365*time_Factor;
  I_0 = 10; M_0 = POPULATION - 0.1 * POPULATION - I_0;
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

  Yp = imatrix(0,No_of_Points, 0, 2);
  Time = dvector(0, No_of_Points);
  Time_Max = dvector(0,No_of_Points); Time_Min = dvector(0,No_of_Points);
  Maxim = dvector(0,No_of_Points); Minim = dvector(0, No_of_Points);
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
  for (i=0;i<Realizations;i++){
    
    SIR_StocRea(Village, &Par, 
		Simlength, modul, 0. , S_I, Time, Yp, I_0);

    save_a_TimeMatrix_to_File_int("sto_", i, Time, Yp, No_of_Points, 3);   
  }
  /* End of STOCHASTIC REALIZATIONS for one particular time */  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  free_imatrix(Yp, 0,No_of_Points, 0, 2);
  free_dvector(Time, 0, No_of_Points);
  free_dvector(Time_Max, 0,No_of_Points); free_dvector(Time_Min, 0,No_of_Points);
  free_dvector(Maxim, 0,No_of_Points); free_dvector(Minim, 0,No_of_Points);

  printf("Fixed Points and Stability of the Closed System (Imm = 0)\n");
  printf("without seasonal forcing\n");
  Fixed_Points(&Par, &f_0, &f_1); Stability(&Par);
  
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}





