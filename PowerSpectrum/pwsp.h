/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	                                                                     */
/*                             David Alonso, 2006 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* XGME_graphs: Parameter to control the graphical execution of the program */
/* #define XGME_graphs */ 
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <nr.h>
#include <nrutil.h>
#include <time.h>
#include <assert.h>
#include <stat.h>

/* Main Functions in main_PW.c */
void ArgumentControl(int, char **);

/* In report_SIR.c */
void modelReport(char *File);

/* In IO.c */ 
void Saving_to_File_powerSpc(char Prefix[], float Time[], float Number[], 
			     int Points, int File_Number);
void Saving_to_File_autocorr(char Prefix[], float Number[], 
			     int Points, int File_Number);
void Showing_powerSpc(float [], float [], int);

/* Functions is powerSpectrum.c */
double total_Power(float *px, int m);
void powerSpectrum_Autocorr(float *x, float *y, int nP);
void powerSpectrum_K_Sets(float *x, float *y, int nP, int k, int m);
void powerSpectrum(float *x, float *y, int nP);
void powerSpectrum_ONE(float *x, int nP);
void accummulating_PowerSpectrum(float *x, float *y, float *px,float *py, int n);
void averaged_PowerSpectrum(double T, float *px, float *py, int n, int Re, int Norma);
void averaged_ONE_PowerSpectrum(double T, float *px, int n, int Re, int Norma);
void accummulating_ONE_PowerSpectrum(float *x, float *px, int n);
void averaged_ONE_Spectral_Density(double T, float *px, int n, int Re, int Norma);
void estimated_Spectral_Density(float *px, float *fr, int n);
void spectral_Density(float data[], float p[], int m, int k, int ovrlap);

/* Function in test_Functions.c */
void sinusoidal_Data(float *X, int, float, double, double, double noise);

/* Functions in autoCorr.c */
void autocorrel(float *p, int n);

/* Functions in autoVar.c */
void autocoVar(float *p, int n);

/* Functions in AR_Process.c */
void AR_p(float *X, int n, int p, float sigma_2, float *alpha);
void AR_p_PWSP_Analytical(float *, int, int, float, float *, float fraction);







