/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                      Computing the Average Power Spectrum                 */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../SIR.h"

int No_of_Points;
int nP; /* Number of points entering the Spectral Analysis */
/* Global Shared parameters main Program <---> ArgumentControl() */
int TIMES;     /* Number of times to be analyzed */
float STEP_SIZE;
float EPSILON;
int Realizations;
int Simlength;
int I_0;       /* Inital number of infective Individuals */
int M_0;      /* Inital number of infective Individuals */
float t_0, t_1;
double a,d,g,mu,b,Imm;
int POPULATION;
double Per;   /* Period of the seasonal Forcing: usually 365 days */
double b_0, b_1;
int TIME_DEPENDENT_PARAM;
int EQ_BETA;

int main(int argc, char **argv)
{
  int N; /* Total Metapopulation Capacity */
  int i,j,k,modul, count_Realizations, spNumber;
  int S_I[3];
  int **Yp;
  double *Time;
  double time_nP;
  float *x, *y, *time;
  float *fr, *px, *py;
  float *fft1, *fft2, *FFT1, *FFT2;
  float T;
  double Psi,Fi; /* Populations fractions at equilibrium from the deterministic system */
  double Factor;
  double Overall_Power, Mean_SQ_Amplitude_x, Mean_SQ_Amplitude_y;
  double Mean_Amplitude_x, Mean_Amplitude_y;
  /* Temporal evolution of the number of strains in each stochastic realization */
  Community Village[No_of_Villages];
  ParamSet Par;

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  Realizations = 5; No_of_Points = 2500; nP = 2048;
  Simlength = 1000000;
  POPULATION = 100000;
  TIMES = 60; STEP_SIZE = 0.25; EPSILON = 0.001;
  InitSeed();
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0., Imm = 1.e-5;
  t_0 = 500.; t_1 = 60000.;
  b_0 = 1.175; b_1 = 0.25; Per = 365.;
  I_0 = 138; M_0 = 44970;
  TIME_DEPENDENT_PARAM = 2; EQ_BETA = 1;
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */

  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);   assert(No_of_Points > nP);

  Yp = imatrix(0,No_of_Points, 0, 2);
  Time = dvector(0, No_of_Points);
  N = POPULATION; /* Village Size */
  settingParameterStruct(&Par);
  modelReport("report.txt");
  /* End of Command line arguments */

  if(TIME_DEPENDENT_PARAM == 2){
    if(EQ_BETA == 1){
      re_setting_Equivalent_BETA(&Par, b_0, b_1);
      TIME_DEPENDENT_PARAM = 0;
      printf("Initial Seasonal Transmission Rate: Beta = %g *(1. + %g)^(Term),\n",
	     b_0, b_1);
      printf("being Term = -1.,1.\n");
      printf("Equivalent Constant Transmission Rate: Beta = %g\tb_1 = 0.\n", Par.Beta);
      //Press_Key();
    }
  }

  if(Simlength < No_of_Points){
    printf("Simulation length must be at least equal to the number of points to be saved\n");
    printf("Simulation length: %d\tPoints to be saved: %d\n", Simlength, No_of_Points);
    exit(0);
  }
  if(Simlength%No_of_Points == 0)
    modul = Simlength/No_of_Points;
  else
    modul = Simlength/No_of_Points + 1;

  x = vector(1,No_of_Points); y = vector(1,No_of_Points);
  time = vector(1, No_of_Points); time_nP = 0.;
  px = vector(0,nP); set_to_value_float(px, nP, 0.);
  py = vector(0,nP); set_to_value_float(py, nP, 0.);
  fft1 = vector(1,2*nP); fft2 = vector(1, 2*nP);
  FFT1 = vector(1,2*nP); set_to_value_float(FFT1, 2*nP, 0.);
  FFT2 = vector(1, 2*nP); set_to_value_float(FFT2, 2*nP, 0.);

  /* Initial Conditions: # of Infective, I_0, and # of Recovered (immune), M_0 */
  printf("Ininitial Condition: Population Values at the stationary state (Imm = 0):\n");
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
  /* Intiating STOCHASTIC REALIZATIONS */
  count_Realizations = 0;
  for (i=0;i<Realizations;i++){
    SIR_StocRea(Village, &Par,
		Simlength, modul, 0., S_I, Time, Yp, I_0);

    if(Yp[No_of_Points-1][1]>0 && Yp[No_of_Points-2][1] > 0){
	count_Realizations++;
	/* Standarizing fluctuations */
	//timeSeries(Time, Yp, No_of_Points, POPULATION, time, x, y, nP, t_0, t_1);
	timeSeries_vanKampen(Psi, Fi, Time, Yp, No_of_Points, POPULATION,
			     time, x, y, nP, t_0, t_1);

	Mean_SQ_Amplitude_x += mean_squared_amplitude(x,No_of_Points);
	Mean_SQ_Amplitude_y += mean_squared_amplitude(y,No_of_Points);

	Mean_Amplitude_x += mean_amplitude(x,nP);
	Mean_Amplitude_y += mean_amplitude(y,nP);

	//Saving_to_File_float_float("x_y", x, y, nP, count_Realizations);
	printf("Standarized fluctuations calculated...\n");

	//powerSpectrum_Autocorr(x, y, nP);                          /* (1) */
	powerSpectrum(x,y,nP);
	accummulating_PowerSpectrum(x,y,px,py,nP/2);                /* (1) */

	printf("Power Spectra calculated and accumulated...\n");
	/* Accumulating last time values... */ time_nP += time[nP];
    }
  }
  /* End of STOCHASTIC REALIZATIONS for one particular time */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  printf("Number of valid Realizations... %d\n", count_Realizations);
  /* Effective average time step... */
  if(count_Realizations > 0){
      T = (time_nP/(double)count_Realizations - time[1])/(float)nP;
      printf("Effective time Sampling interval (days)... %g\n", T);
      Factor = (float)nP*(float)nP*(float)count_Realizations;
      averaged_PowerSpectrum(Factor, px,py, nP/2, count_Realizations, 0); /* (1) */
      /* Here, if last entry is 0, a non-normalized power spectra is calculated */
      Mean_SQ_Amplitude_x /= (double)count_Realizations;
      Mean_SQ_Amplitude_y /= (double)count_Realizations;
      Mean_Amplitude_x /= (double)count_Realizations;
      Mean_Amplitude_y /= (double)count_Realizations;
  }
  /* Defining the analyzed frequencies
     and saving non-standarized power spectra (total power per frequency bin) */
  fr = vector(0,nP);
  for (i=0; i<=nP/2; i++) fr[i] = (float)i/(float)nP/T;
  Saving_to_File_float_float("pwsp", fr, px, nP/2, 0);
  Saving_to_File_float_float("pwsp", fr, py, nP/2, 1);

  Overall_Power = total_Power(px,nP/2);
  printf("Overall_Power = %g\tMean Squared Amplitude = %g\n", Overall_Power, Mean_SQ_Amplitude_x);
  printf("Px(0) = %g\t(Mean_Amplitude_x)^2 = %g\n", px[0], Mean_Amplitude_x*Mean_Amplitude_x);
  Overall_Power = total_Power(py,nP/2);
  printf("Overall_Power = %g\tMean Squared Amplitude = %g\n", Overall_Power, Mean_SQ_Amplitude_y);
  printf("Py(0) = %g\t(Mean_Amplitude_y)^2 = %g\n", py[0], Mean_Amplitude_y*Mean_Amplitude_y);

  /* Calculating and saving Spectral densities (Density Power per unit frequency) */
  estimated_Spectral_Density(px, fr, nP/2);
  Saving_to_File_float_float("spden", fr, px, nP/2, 0);
  estimated_Spectral_Density(py, fr, nP/2);
  Saving_to_File_float_float("spden", fr, py, nP/2, 1);

  Fixed_Points_General(&Par, &Psi, &Fi);
  modelReport("report.txt");

  free_vector(x,1,No_of_Points); free_vector(y,1,No_of_Points);
  free_vector(time,1,No_of_Points);
  free_vector(fr, 0,nP); free_vector(px,0,nP); free_vector(py, 0,nP);
  free_vector(fft1, 1,2*nP); free_vector(fft2, 1,2*nP);
  free_vector(FFT1, 1,2*nP); free_vector(FFT2, 1,2*nP);
  free_imatrix(Yp, 0,No_of_Points, 0, 2);
  free_dvector(Time, 0, No_of_Points);

  printf("\nEnd of progam\n");
  return (0);
}
