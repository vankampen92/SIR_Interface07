/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "pwsp.h"

double total_Power(float *px, int m)
{
  int i;
  double power;

  power=0.;
  for(i=0; i<=m; i++) power += px[i];

  return power;
}

void powerSpectrum_ONE(float *x, int nP)
{
  /* Power Spectra of one series x */
  int i,j,k;
  unsigned long n;
  float factor;
  double *fft1;

  fft1 = dvector(1,2*nP);
  n = nP;

  i=1;
  for(k=1; k<=2*nP; k++){
    if(k%2 == 1){
         fft1[k] = x[i];
	 i++;
      }
      else{
	 fft1[k] = 0.;
      }
  }

  dfour1(fft1, nP, 1);

  j=0;
  for(i=1; i<2*n; i+=2){
    j++;
    x[j] = fft1[i]*fft1[i] + fft1[i+1]*fft1[i+1];
  }

  for(i=1; i<n/2; i++){
    x[i+1] = (x[i+1] + x[n-i+1]);
  }

  free_dvector(fft1, 1,2*nP);
}

void averaged_ONE_PowerSpectrum(double Factor, float *px, int n, int Re, int Norma)
{
    /* The estimated discret FFT is related to the its continuous integral
       counterpart by a factor of T, the dicretization time interval */
  int i;

  /* Factor = (float)Re * 4.*(float)n *(float)n / (float)T; */
  /* This is a normalizing  factor */
  for(i=0; i<=n; i++){
    px[i] /= Factor;
  }

  /* Normalizing */
  if(Norma == 1){
    Norma_Total_a_J(px, 0, n-1);
  }
}

void averaged_ONE_Spectral_Density(double Factor, float *px, int n, int Re, int Norma)
{
    /* The estimated discret FFT is related to the its continuous integral
       counterpart by a factor of T, the dicretization time interval */
  int i;

  /* This is a normalizing  factor */
  for(i=0; i<=n; i++){
    px[i] /= Factor;
  }

  /* Normalizing */
  if(Norma == 1){
    Norma_Total_a_J(px, 0, n-1);
  }
}

void estimated_Spectral_Density(float *px, float *fr, int n)
{
  /* px[i] are not densities, but the total power associated to the bin
     fr[i]-fr[i-1]. This routine will turn px[] into a truely spectral
     density
  */
  int i;

  for(i=0; i<n; i++)
    px[i] = px[i]/(fr[i+1]-fr[i]);
}

void accummulating_ONE_PowerSpectrum(float *x, float *px, int n)
{
  int i;

  for(i=0; i<=n; i++)  px[i] += x[i+1];
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */
/* Simultaneous calculation of the power spectra from two time series */
void powerSpectrum_K_Sets(float *x, float *y, int nP, int k, int m)
{
  /* . k is the number of subsets used to take an average power spectral
     density.
     . nP MUST be an integer power of two
  */
  float *px, *py;
  int i;

  if(k == 1){

    powerSpectrum(x,y,nP);

  }
  else{

    if(nP != 4*k*m){
      printf("The number of temporal points points is not correct!!\n");
      printf("Time Series Length is %d, but should be %d\n", nP, 4*k*m);
      exit(0);
    }

    px = vector(1,m+1); py = vector(1,m+1);

    spectral_Density(x, px, m, k, 0);
    spectral_Density(y, py, m, k, 0);

    for(i=1; i<=m; i++) {
      x[i] = px[i];
      y[i] = py[i];
    }

    free_vector(px, 1,m+1); free_vector(py, 1,m+1);
  }
}

void powerSpectrum(float *x, float *y, int nP)
{
  /* Power Spectra of two series x and y at at time */
  /* nP must be an interger power of 2 */

  int i,j;
  unsigned long n;
  float factor;
  float *fft1, *fft2;

  fft1 = vector(1,2*nP); fft2 = vector(1, 2*nP);
  n = nP;

  twofft(x,y,fft1,fft2,n);

  j=0;
  for(i=1; i<2*n; i+=2){
    j++;
    x[j] = fft1[i]*fft1[i] + fft1[i+1]*fft1[i+1];
    y[j] = fft2[i]*fft2[i] + fft2[i+1]*fft2[i+1];
  }

  for(i=1; i<n/2; i++){
    x[i+1] = x[i+1] + x[n-i+1];
    y[i+1] = y[i+1] + y[n-i+1];
  }

  free_vector(fft1, 1,2*nP); free_vector(fft2, 1,2*nP);
}

void accummulating_PowerSpectrum(float *x, float *y,
				 float *px,float *py, int n)
{
  int i;
  for(i=0; i<=n; i++){
    px[i] += x[i+1];
    py[i] += y[i+1];
  }
}

void averaged_PowerSpectrum(double Factor, float *px, float *py, int n, int Re, int Norma)
{
    /* The estimated discret FFT is related to the its continuous integral
       counterpart by a factor of T, the dicretization time interval */
  int i;
  /* Factor is defined in the  main program */
  // Factor = 4.*(float)n*(float)n * (float)Re / (float)T;

  for(i=0; i<=n; i++){
    px[i] /= Factor;
    py[i] /= Factor;
  }

  /* Normalizing */
  if(Norma == 1){
    Norma_Total_a_J(px, 0, n-1);
    Norma_Total_a_J(py, 0, n-1);
  }
}

void powerSpectrum_Autocorr(float *x, float *y, int nP)
{
  /* The Power Spectrum is calculated as the Fourier trasnform of the
     autocorrelation function */

  int i,k;
  unsigned long n;
  double *fft1, *fft2;
  n = nP;

  fft1 = dvector(1, 2*nP); fft2 = dvector(1, 2*nP);

  /* Calculating first the correlations functions of x and y series */
  autocoVar(x, nP); autocoVar(y, nP);

  i=1;
  for(k=1; k<=2*nP; k++){
    if(k%2 == 1){
         fft1[k] = x[i];
         fft2[k] = y[i];
	 i++;
      }
      else{
	 fft1[k] = 0.;
	 fft2[k] = 0.;
      }
  }
  dfour1(fft1, nP, 1);  dfour1(fft2, nP, 1);

  for(i=0; i<nP/2; i++){
      x[i+1] = fft1[1+2*i];
      y[i+1] = fft2[1+2*i];
  }

  free_dvector(fft1, 1, 2*nP); free_dvector(fft2, 1, 2*nP);
}
