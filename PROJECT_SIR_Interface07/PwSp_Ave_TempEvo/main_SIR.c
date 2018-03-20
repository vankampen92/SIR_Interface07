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
/* Global Shared parameters main Program <---> ArgumentControl() */
int TIMES;     /* Number of times to be analyzed */
int TRANSIENT;
int No_of_Points;
int k_Sets;
float STEP_SIZE; 
float EPSILON;
int Realizations;
int Simlength;
int I_0, M_0;      /* Inital number of infective Individuals */
double a,d,g,mu,b,Imm;
int POPULATION;
float t_0, t_1;
double Per;  /* Period of the seasonal Forcing: usually 365 days */
double b_0, b_1;
int TIME_DEPENDENT_PARAM;
int DISCARDING_EXTINCTIONS;
int EQ_BETA;
float time_Factor;
void fill_FileNames(char Pref[2][25], int);

int main(int argc, char **argv)
{
  int N;          /* Total Metapopulation Capacity */
  int mm, no_Fr; /* Number of frequency points to be saved */
  int i,j,jj,k,modul, spNumber, Bad_Times, Extinction; 
  int Yp[1][3], Strain[3]; /*0: S; 1: I; 2: R */
  double Time[1];
  double *Time_values;
  float *x, *y, *time;
  float *fr, *px, *py;
  int count_Realizations;
  /* Time_values[TIMES] = 
     {0.2, 0.3, 0.4, 0.5, 0.7, 1., 2., 5., 6., 8., 10., 12., 15., 20.}; */
  /* Temporal evolution of the number of strains in each stochastic realization */
  int new; /* New infections during the time interval */
  char Files[2][25];
  FILE *FP; char file[12]; FILE *Fp; char File[12];
  double Psi,Fi; /* Populations fractions at equilibrium from the deterministic system */ 
  double Factor;
  double Overall_Power, Mean_SQ_Amplitude_x, Mean_SQ_Amplitude_y;
  double Mean_Amplitude_x, Mean_Amplitude_y;
  Community Village[No_of_Villages];
  ParamSet Par;

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  TIME_DEPENDENT_PARAM = 0; EQ_BETA = 1; DISCARDING_EXTINCTIONS = 0;
  Realizations = 1; No_of_Points =  1000;
  Simlength = 1000; time_Factor = 1.;
  POPULATION = 10000;
  TIMES = 1200; TRANSIENT = 100; STEP_SIZE = 14.; EPSILON = 0.1; k_Sets = 1;
  InitSeed();
  b = 1.175; a = 0.; d = 0.000055; g = 1./13.; mu = 0.; Imm = 0.00001;
  b_0 = 1.175; b_1 = 0.25; Per = 365.;
  I_0 = 10; M_0 = 0.1 * POPULATION;
  t_0 = 0.; t_1 = 0.; 
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */ 
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);
  
  mm = TIMES/4/k_Sets;
  if(k_Sets > 1) no_Fr = mm;
  else           no_Fr = TIMES/2;          
  printf("The number of frequency points will be %d\n", no_Fr);
  
  Time_values = dvector(0,TIMES+TRANSIENT-1);
  N = POPULATION; 
  settingParameterStruct(&Par);
  modelReport("report.txt");
  modul = Simlength/No_of_Points;
  
  if(t_1 > 0.)
    STEP_SIZE = (t_1-t_0)/(float)(TIMES+TRANSIENT);
  else
    t_1 = t_0 + STEP_SIZE * (TIMES+TRANSIENT);
  
  if(TIME_DEPENDENT_PARAM == 2){
    if(EQ_BETA == 1){
      re_setting_Equivalent_BETA(&Par, b_0, b_1);
      TIME_DEPENDENT_PARAM = 0;
      printf("\n  Non Seasonal Forcing with an equivalent effective Transmission Rate (Beta)\n"); 
      printf("  Initial Seasonal Transmission Rate: Beta = %g *(1. + %g)^(Term),\n", 
	     b_0, b_1);
      printf("  being Term = -1.,1.\n"); 
      printf("  Equivalent Constant Transmission Rate: Beta = %g\tb_1 = 0.\n", Par.Beta);
      //Press_Key();
    }   
  }

  printf("Sampling Interval (days): %g\n", STEP_SIZE);
  printf("Entering Generation of Stochastic Realizations...\n");
  //Press_Key();
  /*setting times*/
  for(i=0; i<TIMES+TRANSIENT; i++)
    Time_values[i] = t_0 + (double)(i+1) * STEP_SIZE;
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */  
  /* Main loop: For any final time, a number of REALIZATIONS is computed */
  
  x = vector(1,TIMES); y = vector(1,TIMES); 
  time = vector(1, TIMES); 
  px = vector(0,TIMES); set_to_value_float(px, TIMES, 0.);  
  py = vector(0,TIMES); set_to_value_float(py, TIMES, 0.); 

  /* Initial Conditions: # of Infective, I_0, and # of Recovered (immune), M_0 
     Notice that these are for an equivalent Beta  */
  printf("Ininitial Condition: Population Values at the stationary state (Imm = 0):\n");
  re_setting_Equivalent_BETA(&Par, b_0, b_1);
  
  Fixed_Points(&Par, &Psi, &Fi);  
  I_0 = POPULATION*Fi; M_0 = POPULATION*(1.-Fi-Psi);
  printf("Susceptible, S(0) = %f\tInfective, I(0) = %d\n", POPULATION*Psi, I_0); 
  printf("\n");
  printf("Initial Condition: Population Values at the statinary state (Imm = %g):\n",
	 Par.Imm);
  Fixed_Points_General(&Par, &Psi, &Fi);
  I_0 = POPULATION*Fi; M_0 = POPULATION*(1.-Fi-Psi);
  printf("Susceptible, S(0) = %f\tInfective, I(0) = %d\n", POPULATION*Psi, I_0); 
  //Press_Key();
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  Mean_SQ_Amplitude_x = 0; Mean_SQ_Amplitude_y = 0;
  Mean_Amplitude_x = 0.;   Mean_Amplitude_y = 0;
  count_Realizations = 0;
  /* Initiating Stochastic Realizations */
  for (i=0;i<Realizations;i++){
    Time[0] = t_0;
    Bad_Times = 0; Extinction = 0;
    jj = 0.; j = 0;
    while (j < (TIMES+TRANSIENT) && Extinction == 0 && Bad_Times == 0){

      SIR(Village, &Par, 
	  Simlength, modul, Time_values[j], Strain, Time, Yp, &new);
      
      if(DISCARDING_EXTINCTIONS == 1){
	if(Strain[1] == 0) Extinction++;
      }
      if((Time[0] > Time_values[j] - EPSILON) && Extinction == 0){
	
	if(j >= TRANSIENT){
	  jj++;
	  x[jj] = ((float)Strain[0]-(float)POPULATION*(float)Psi)/sqrt((double)POPULATION);
	  y[jj] = ((float)Strain[1]-(float)POPULATION*(float)Fi)/sqrt((double)POPULATION);
	  time[jj] = Time[0];
	  /* Time[0] is always the last time which is the closest but lower to Time_values[j] */
	  
	}
      }
      else{
	Bad_Times++;
      }
      j++;
      /* go further to the next time */
    }

    printf("There were %d extinctions\tThese series are discarded to calculate power spectra\n", 
	   Extinction);
    printf("Time failed if Bad_Times is 1\t Bad_Times = %d\n", Bad_Times);
    printf("   If you don't like to fail, try to choose a larger EPSILON -E %g,\n",EPSILON);
    printf("   because accuracy might be too stringent (EPSILON too low)!\n");
    
    /* Calculating Spectra: 
       Times series without missing values */
    if(Bad_Times == 0 && Extinction == 0){
      Mean_SQ_Amplitude_x += mean_squared_amplitude(x,TIMES); 
      Mean_SQ_Amplitude_y += mean_squared_amplitude(y,TIMES);
	
      Mean_Amplitude_x += mean_amplitude(x,TIMES);
      Mean_Amplitude_y += mean_amplitude(y,TIMES);
      
      powerSpectrum_K_Sets(x, y, TIMES, k_Sets, mm);         
      accummulating_PowerSpectrum(x,y,px,py,no_Fr);                /* (1) */

      count_Realizations++;
    }
    printf("Realization: %d of a total of %d\n", i+1, Realizations);
  }
  
  /* End of STOCHASTIC REALIZATIONS */
  /* Saving Power Spectra */
  printf("Power Spectra calculated over %d stochastic realizations\n", count_Realizations);
  if(count_Realizations > 0){
    if(k_Sets == 1)
       Factor = 4.*(float)no_Fr*(float)no_Fr*(float)count_Realizations;
    else
       Factor = (float)count_Realizations;
 
    averaged_PowerSpectrum(Factor, px,py, no_Fr, count_Realizations, 0); /* (1) */

    printf("Power Spectra calculated and accumulated over %d stochastic realizations\n", 
	   count_Realizations);
    Mean_SQ_Amplitude_x /= (double)count_Realizations;
    Mean_SQ_Amplitude_y /= (double)count_Realizations; 
    Mean_Amplitude_x /= (double)count_Realizations;
    Mean_Amplitude_y /= (double)count_Realizations;
  }
  /* Defining the analyzed frequencies
     and saving non-standarized power spectra (total power per frequency bin) */
  fr = vector(0,TIMES);
  //for (i=0; i<=no_Fr; i++) fr[i] = (float)i/2./(float)no_Fr/STEP_SIZE;
  for (i=0; i<=no_Fr; i++) fr[i] = (float)i/2./(float)no_Fr/STEP_SIZE * 365.;
  
  Saving_to_File_float_float("pwsp", fr, px, no_Fr, 0, 0);
  Saving_to_File_float_float("pwsp", fr, py, no_Fr, 1, 0);

  Overall_Power = total_Power(px,no_Fr);
  
  printf("Overall_Power = %g\tMean Squared Amplitude = %g\n", Overall_Power, Mean_SQ_Amplitude_x);
  printf("Px(0) = %g\t(Mean_Amplitude_x)^2 = %g\n", px[0], Mean_Amplitude_x*Mean_Amplitude_x);
  Overall_Power = total_Power(py,no_Fr);
  printf("Overall_Power = %g\tMean Squared Amplitude = %g\n", Overall_Power, Mean_SQ_Amplitude_y);
  printf("Py(0) = %g\t(Mean_Amplitude_y)^2 = %g\n", py[0], Mean_Amplitude_y*Mean_Amplitude_y);
  /* Calculating and saving Spectral densities (Density Power per unit frequency) 
     in a way that the sum over the frequency grid produces the Mean Squared Amplitude */
  estimated_Spectral_Density(px, fr, no_Fr);  
  Saving_to_File_float_float("spden", fr, px, no_Fr, 0, 0);

  estimated_Spectral_Density(py, fr, no_Fr);
  Saving_to_File_float_float("spden", fr, py, no_Fr, 1, 0);
  
  Fixed_Points_General(&Par, &Psi, &Fi);  

  modelReport("report.txt");
  free_vector(x,1,TIMES); free_vector(y,1,TIMES); 
  free_vector(time,1,TIMES);
  free_vector(fr, 0,TIMES); free_vector(px,0,TIMES); free_vector(py, 0,TIMES);
  free_dvector(Time_values, 0,TIMES+TRANSIENT-1);
  
  printf("\nEnd of progam\n");
  return (0);
}







