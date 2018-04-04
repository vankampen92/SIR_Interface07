/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	       (implementing spctrm( ... )form Numerical Recipes)            */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../pwsp.h"

/* Global Shared parameters main Program <---> ArgumentControl() */
int No_of_Points;
int REPLICAS;
int ONE;

int main(int argc, char **argv)
     /* 
	This main program calculates a spectral estimate of a data set by reading 
	it in as segments, windowning, Fourier transforming, and accumulating the 
	power spectrum. Data segments may or may not be overlapped at the decision of
	the user. The test data set is a AR(p) process. The analytical predicted
	power spectrum is also calculated for comparison purposes.
     */	
{
  float fraction;
  double Factor;
  double Overall_Power, Mean_SQ_Amplitude;
  FILE *fp;
  float *fr;
  float *fi;
  float *X, *px;
  float *Y, *py;
  int i, k,j, p;
  float sigma;

  /* Default parameter values */
  p = 2; ONE = 0; 
  No_of_Points = 256;
  REPLICAS = 1000;

  if(argc>1) 
    ArgumentControl(argc,argv); 

  fi = vector(1,p);
  /* Initializing parameters */
  //fi[1] = 1.318; fi[2] = -0.1634; sigma = 289.2; /* Sun spost */
  fi[1] = 1.8645; fi[2] = -0.9608; sigma = 0.25;
  fr = vector(0,No_of_Points);
  X = vector(1,No_of_Points);  px = vector(0, No_of_Points); 
  Y = vector(1,No_of_Points);  py = vector(0, No_of_Points); 
  set_to_value_float(px, No_of_Points, 0.);  
  set_to_value_float(py, No_of_Points, 0.);

  /* Averaging over REPLICAS */
  printf("Realization... \n");
  for(i=0; i<REPLICAS; i++){
    /* AR(p): Data Generation */
    AR_p(X, No_of_Points, p, sigma, fi);  
    AR_p(Y, No_of_Points, p, sigma, fi);
    if(ONE == 1){
      powerSpectrum_ONE(X, No_of_Points);
      accummulating_ONE_PowerSpectrum(X, px, No_of_Points/2);
    }
    else{
      powerSpectrum(X,Y,No_of_Points);
      accummulating_PowerSpectrum(X,Y,px,py,No_of_Points/2); /* (1) */
    }
    printf("%d...", i);
  }
  /* End Realizations */
  for (i=0; i<=No_of_Points; i++) fr[i] = (float)i/(float)No_of_Points;
  
  Factor = (float)No_of_Points*(float)No_of_Points *(float)REPLICAS; 
  if(ONE == 1){
    averaged_ONE_PowerSpectrum(Factor, px, No_of_Points/2, REPLICAS, 0);
    Saving_to_File_float_float("pwsp", fr, px, No_of_Points/2, 2);
  }
  else{
    averaged_PowerSpectrum(Factor, px,py, No_of_Points/2, REPLICAS, 0); /* (1) */
    Saving_to_File_float_float("pwsp", fr, px, No_of_Points/2, 2);
    Saving_to_File_float_float("pwsp", fr, py, No_of_Points/2, 4);
  }
  printf("Done!!!\n");
  /* End of Averaging over REPLICAS */

  AR_p(X, No_of_Points, p, sigma, fi);  
  Mean_SQ_Amplitude = mean_squared_amplitude(X,No_of_Points);
  Overall_Power = total_Power(px,No_of_Points/2);
  printf("Overall_Power = %g\tMean Squared Amplitude = %g\n", Overall_Power, Mean_SQ_Amplitude);

  estimated_Spectral_Density(px, fr, No_of_Points/2);
  Saving_to_File_float_float("spden", fr, px, No_of_Points/2, 4);

  /* Analytical Predicted Power Spectrum for an AR(2) */
  fraction = 0.5;
  for (i=0; i<=No_of_Points; i++) fr[i] = (float)i*fraction/(float)No_of_Points;
  AR_p_PWSP_Analytical(px, No_of_Points, p, sigma, fi, fraction);
  Saving_to_File_float_float("pwAn", fr, px, No_of_Points/2, 0);

  free_vector(fr, 0, No_of_Points); free_vector(fi, 1,p);
  free_vector(px, 0, No_of_Points); free_vector(X, 1, No_of_Points); 
  free_vector(py, 0, No_of_Points); free_vector(Y, 1, No_of_Points);
  printf("All power spectra has been estimated successfully!!\n");
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}





