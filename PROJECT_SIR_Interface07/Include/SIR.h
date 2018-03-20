/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <nr.h>
#include <nrutil.h>
#include <time.h>
#include <assert.h>
#include <stat.h> 

#define EVOLUTION_PARTICULAR_PATCHES 0 /* 1/0. When 1 the patch content is shown 
					  and saved for every particular patch 
				       */
typedef struct communityinfo
{
  int N;             /* Total Patch capacity: Human Constant Population     */
  int m;             /* Total human population of susceptibles m = I[0]     */
  int n;             /* n = I[1]  Number of infected                        */
  int * I;          
  double ratePatch;  /* Transition probability at this patch                    */
  double * rate;     /* Transition probability for an individual at each class  */
  double * rToI; 
  double H_infect;
  struct point position;  /* Patch Geographical Coordinates  */ 
}Community;
  
typedef struct parameterinfo
{
  double Beta;
  double Alpha;
  double Gamma;
  double Delta;
  double Mu;
  double Imm;
}ParamSet;

typedef struct Parameter_Tableinfo
{
  int POPULATION;    /* System Max Capicity: Sum of all max population per patch */

  int TIME_DEPENDENCE; 
  int DISCARDING_EXTINCTIONS;

  ParamSet * Model_Parameters;
  int no_Patch;        
  Community ** Patch;

  int Simulation_Length; 

  double Time_0;         
  int Infectious_Time_0;     
  int Susceptible_Time_0;    
  double Time_Scale_Unit;
    
  double Time;                     /* Current time of the system                           */             
  double ** Populations_Per_Patch; /* Population types per path and topal population types 
				      in the system at a given time, Time                  */
  int * Total_Populations;        /* Total population per type at a given time, Time      */   
}Parameter_Table;

/* In Main Functions defining the simulation models and cotroling argument Input */
int SIR( Parameter_Table * Table, double Time_Initial, double * Time_Final, int * New );

void SIR_StocRea( Community *Village, ParamSet *Set, 
		  int Simlength, int modul, double t_final, 
		  int Sum_I[], double *Time, int **, int I_0 );

void ArgumentControl(int, char **);

void poincareMap(int i, int j, float Pop[], int N, float x_i,float x_s, 
		 int No_of_Points, void (*)(double,double[],double[]),
		 double *Time, double **X);

void integration_Double(float [], int S, float,float, int, void (*)(double,double[],double[]),
		    double *, double **);

void integration(float [], int,float,float,int, void (*)(float,float[],float[]));
void deriva (double, double[], double[]);
void derivaSIR (double, double[], double[]);
void derivaSEIR (double, double[], double[]);

/* In init_SIR.c */
void Init_SumPop(Community *, int [], int);
void settingParameterStruct(ParamSet *Par);
void Initial_Condition(Community *, int, int);
void changingTimeScale(ParamSet *Par, double TimeScale);
void re_settingParamStruct(ParamSet *Par, double d_Rel, double b_Rel);
void re_setting_Equivalent_BETA(ParamSet *Par, double b_0, double b_1);

/* In report_SIR.c */
void modelReport(char *File);

/* Functions in func_SIR.c */
int which_Village(double max_Probability, Community *Pop);
struct point RandomPosition(int Length);
int Total_Number_of_Infective(int *SumPop, int Sp);
void Fixed_Points(ParamSet *P, double *Psi, double *Fi);
void Stability(ParamSet *P);
float approximate_Endogenous_Period(ParamSet *P);

/* Functions in IO_SIR.c */
void ShowVillageContent(Community *pVil, int i);
void ShowVillage_StrainEvolution(Community *pVil, int i, double time);
void SaveVillageContent(char Pre[], Community *pVil, int, int, double);
void SaveTimeSeries(Community *pVil, int i, double te);
void write_vectors(float *time, float *x, float *y, int nP);

/* Functions in even_TimeSeries.c */
void pseudoEven_timeSeries(double *Time, int **Yp, int NoP, int nPop, 
			   float *time, float *x, float *y, int n, 
			   float t_0, float t_1);
void Even_timeSeries(double *Time, int **Yp, int NoP, int nPop, 
		     float *time, float *x, float *y, int n, 
		     float t_0, float t_1);
void timeSeries(double *Time, int **Yp, int NoP, int nPop, 
		float *time, float *x, float *y, int n, 
		float t_0, float t_1);

/* Functions in model_SIR.c */
int execute_Step(Community *Pop, ParamSet *Se, double *Total_Rate, 
		    double Total_Probability, int SumPop[]);
void Temporal_Dynamics(Community *Village, ParamSet *Se, 
		       double *Rate, double *Total_Probability, int SumPop[]);

/* Functions in evolving_Parameters */
void evolving_Parameters(double time, float Factor, ParamSet *p, int Typus);
double Time_of_the_Year(double time);
double Time_of_the_Year_Two_Terms(double time);

/* Function in sampling_Parameters.c */
void sample_Parameters(ParamSet *P, 
		       float *range_R, float *range_L, float *range_D);

/* Function in peak_To_trough_Ratio.c */
float peak_to_Trough_Global_Range(double *Time, int **Yp, int n);
void peak_to_Trough_Maxims_Minims(double *Time, int **Yp, 
				  int No, int Level,
				  double *Time_Max, double *Time_Min, 
				  double *Maxim, double *Minim, int *noVal);
void saving_peak_to_Trough(double *Time_Max, double *Time_Min, 
			   double *Maxim, double *Minim, int No);
float peak_to_Trough_Average_Ratio(double *Maxim, double *Minim, int No);
float Relative_Amplification_Ratio(double *Maxim, int No, int I_0);

/* Functions in det_seasonal_amplif.c */
float det_season_amplif(int MSD, ParamSet *Pa, float T_sup, int no_of_Maxima, int SI);

/* Functions in sto_seasonal_amplif.c */
void the_Largest_Maxima(double *Maxim, int no, int no_of_Max);
float sto_season_amplif(int MSD, ParamSet *Pa, float T_sup, int no_of_Maxima, int SI);

/* Function in relativeAmplification_forcing.c */
float Relative_Amplification_Ratio_Forcing(double *Time_Max, double *Maxim, int No,
					   int N, double f_0, double f_1,
					   ParamSet *Par, int No_of_Points);

/* Functions is powerSpectrum.c */
double total_Power(float *px, int m);
void powerSpectrum_Autocorr(float *x, float *y, int nP);
void powerSpectrum(float *x, float *y, int nP);
void powerSpectrum_ONE(float *x, int nP);
void accummulating_PowerSpectrum(float *x, float *y, float *px,float *py, int n);
void averaged_PowerSpectrum(double T, float *px, float *py, int n, int Re, int Norma);
void averaged_ONE_PowerSpectrum(double T, float *px, int n, int Re, int Norma);
void accummulating_ONE_PowerSpectrum(float *x, float *px, int n);
void averaged_ONE_Spectral_Density(double T, float *px, int n, int Re, int Norma);
void estimated_Spectral_Density(float *px, float *fr, int n);

/* Functions in timeSeries_VanKampen.c */
void timeSeries_vanKampen(double Ave_x, double Ave_y,
			  double *Time, int **Yp, int NoP, int nPop, 
			  float *time, float *x, float *y, int n, 
			  float t_0, float t_1);

/* Functions in autoCorr.c */
void autocorrel(float *p, int n);

/* Functions in ./Fig_eps/append_Fig_TEX.c */
void append_Fig_TEX(FILE *F_TEX, char *name, float imm);

/* In default_values.c... */
void default_values_Measles_01();
void default_values_Measles_02();
void values__Measles_02_SEIR();
void values__Measles_02_SIR();


