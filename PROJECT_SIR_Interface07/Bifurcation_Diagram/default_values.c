/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                          SIR MODEL: Default Values                        */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                            David Alonso, 2000 (c)                         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../SIR.h"
#include "../SIR_SEIR.h"
extern float x_i, x_s;
extern float p_1, p_2;
extern double a;
extern double mu;
extern double b;
extern double b_0, b_1;
extern double Per;
extern double si;
extern double g;
extern double d;
extern double Imm;
extern int N;     /* Number of equations */
extern int No_of_Points; 
extern int TIME_DEPENDENT_PARAM;
extern int POPULATION;
extern float timeFactor;
/* timeFactor is defined as the 1/U, where U is temporal unit used in days */

extern double b_m, b_M;
extern int No_of_Beta;
extern int No_of_Lattice;
extern int SEIR, sir;
extern int MODEL;

extern ParamSEIR *Pe;
extern ParamSIR *Pi;

void default_values_Measles_01()
{
  /* Default values in days^(-1)*/
  POPULATION = 100000;
  No_of_Points = 1000; TIME_DEPENDENT_PARAM = 0;
  x_i=0.; x_s=36500.;
  p_1=0.1; p_2=0.0001;
  b = 1.175; a = 0.; d = 5.5e-05; g = 1./13.; mu = 0.; Imm = 1.e-05;
  b_0 = b; b_1 = 0.25; timeFactor = 1.; Per = 365.*timeFactor;
}

void default_values_Measles_02()
{
  /* Default values in years^(-1), by default the model is SIR */
  MODEL = 0;
  POPULATION = 1000000;  SEIR = 0; sir = 0;
  TIME_DEPENDENT_PARAM = 3;
  x_i=0.; x_s=500.;
  No_of_Points = (int)x_s; /* One only data per year (the forcing period) in order to calculate 
			      for the Poincare Map */ 
  p_1=0.1; p_2=0.0001;
  b = 365. * 1.175; a = 0.; d = 0.02; g = 365. * 1./13.; mu = 0.; Imm = 365. * 1.e-06;
  si = 365 * 1./8.;
  b_0 = b; b_1 = 0.25; timeFactor = 1./365.; Per = 365.*timeFactor;
  No_of_Lattice = 40;
  No_of_Beta = 100;
  b_m = 200.; b_M = 2000.; 
}

void values_Measles_02_SEIR(ParamSEIR *P)
{
  Pe = P;       /* Pe will point to the parameter table 
		   defined in the main program */
 
  si = 365. * 1./8.; g = 365. * 1./5.;
  P->b = b;     /* Basal Transmision rate */
  P->b_0 = b;   /* Basal Transmision rate */
  P->b_1 = b_1; /* Seasonality Parameter  */
  P->si = si;   /* 1./si, average time in the Exposed, but no infected phase */
  P->g  = g;    /* 1./g, average time in the Infectious phase */
  P->nu = d;    /* Per capita birth rate */
  P->d = d;     /* Per capita mortality rate */
  P->a = a;     /* Disease-induced death rate */
  P->mu = mu;   /* Rate of loss of immunity */
  P->Imm = Imm; /* Immigration parameter */
  P->Per = Per; /* Period of seasonal forcing */
}

void values_Measles_02_SIR(ParamSIR *P)
{
  Pi = P;       /* Pi will point to the parameter table 
		   defined in the main program and for which
		   the there has been space reserved */
  
  P->b = b;     /* Basal Transmision rate */
  P->b_0 = b;   /* Basal Transmision rate */
  P->b_1 = b_1; /* Seasonality Parameter  */
  P->g  = g;    /* 1./g, average time in the Infectious phase */
  P->nu = d;    /* Per capita birth rate */
  P->d = d;     /* Per capita mortality rate */
  P->a = a;     /* Disease-induced death rate */
  P->mu = mu;   /* Rate of loss of immunity */
  P->Imm = Imm; /* Immigration parameter */
  P->Per = Per; /* Period of seasonal forcing */
}






