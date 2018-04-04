/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	               (implementing spectral_Density( ... )                 */ 
/*	            based on spctrm(..) form Numerical Recipes)              */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../pwsp.h"

         /* Global Shared parameters main Program <---> ArgumentControl() */
int p;  /* Order of an AR(p) process                                     */
int m; /* Number of different frequencies in the PS generation...       */
int ovrlap;
char *File;

int main(int argc, char **argv)
     /* 
	This main program calculates a spectral estimate of a data set by reading 
	it in as segments, windowning, Fourier transforming, and accumulating the 
	power spectrum. Data segments may or may not be overlapped at the decision of
	the user. The test data set is a AR(p) process. The analytical predicted
	power spectrum is also calculated for comparison purposes.
     */	
{
  int No_of_Points;
  float fraction;
  FILE *fp;
  float *pe, *fr;
  float *X, *px, *fi;
  int i, k,j;
  float sigma;
  char Fi[100] = "timeData.dat";

  /* Default parameter values */
  File = Fi;   p = 2;   
  ovrlap = 1;  m = 256; k = 10; 

  if(argc>1) 
    ArgumentControl(argc,argv); 

  if(ovrlap == 1)
    No_of_Points = (2*k + 1) * m;
  else 
    No_of_Points = 4*k*m;

  fi = vector(1,p);
  /* Initializing parameters */
  //fi[1] = 1.318; fi[2] = -0.1634; sigma = 289.2; /* Sun spost */
  fi[1] = 1.8645; fi[2] = -0.9608; sigma = 0.25;
  
  pe = vector(1,m+1); px = vector(0,m); fr = vector(0,m);
  X = vector(1,No_of_Points);

  /* AR(p): Data Generation */
  fp = fopen(File, "w");
  AR_p(X, No_of_Points, p, sigma, fi);     /* AR process of simga^2 = 1, and 
					      coeficients all equal to 0.5 */ 
  for(i=1; i<=No_of_Points; i++)
    fprintf(fp, "%f\n", X[i]);
  fclose(fp);
  /* End of Data Generation */

  fp = fopen(File, "r");
  printf("Calculating and saving first power spectrum...\n");
  //spctrm(fp, pe, m, k, ovrlap);
  spectral_Density(X, pe, m, k, ovrlap);

  for(i=0; i<=m; i++) {
    fr[i] = (float)i/(2. * (float)m);
    px[i] = pe[i+1];
  }
  Saving_to_File_float_float("pwsp", fr, px, m, ovrlap); 
  printf("Done!!!\n");
  fclose(fp);

  free_vector(pe, 1,m+1); free_vector(px, 0,m); free_vector(fr, 0, m); 
  free_vector(X, 1, No_of_Points); 

  /* Analytical Predicted Power Spectrum for an AR(2) */
  fraction = 0.5;
  fr = vector(0, No_of_Points); px = vector(0, No_of_Points);
  for (i=0; i<=No_of_Points; i++) fr[i] = (float)i*fraction/(float)No_of_Points;
  AR_p_PWSP_Analytical(px, No_of_Points, p, sigma, fi, fraction);
  Saving_to_File_float_float("pwAn", fr, px, No_of_Points/2, 0);

  free_vector(fr, 0, No_of_Points); 
  free_vector(px, 0, No_of_Points); free_vector(fi, 1,p);

  printf("All power spectra has been estimated successfully!!\n");
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}





