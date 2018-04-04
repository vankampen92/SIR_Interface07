/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*              implementing void period( ... ) from Numerical Recipes.      */
/*                       Uneven distribution of temporal data.                */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../pwsp.h"
/* Global Shared parameters main Program <---> ArgumentControl() */
int No_of_Points;
char *File;
float t_0;

int main(int argc, char **argv)
{
  float *fr; float *p;
  float *x;   float *y; float *time;
  int nP,nF;
  int i;
  //unsigned long nout, jmax; /* When using fasper() */
  int nout, jmax;
  float prob;
  float hifac, ofac, T; /* Time Spanning of the analyzed Time Series */
  float fr_c;
  char Fi[100] = "timeData.dat";

  /* Defauld values */
  File = Fi;  t_0 = 1.;
  ofac = 4.;  No_of_Points = 5000;

  if(argc>1) 
    ArgumentControl(argc,argv); 

  /* Preparing x and y after discarding transients */
  time = vector(1,No_of_Points);
  x = vector(1,No_of_Points); y = vector(1, No_of_Points);
  /* Selecting the stationary part of a 2d time series
     from a (time, x, y) three column file */
  timeSeries(File, time, x, y, No_of_Points, &nP, t_0);
  
  hifac = 0.5; /* Exploring up to hifac*Nyquist_Frequency */
  T = time[nP] - time[1];
  fr_c = (float)nP/2/T;
  nF = 0.5*ofac*hifac*nP;
  printf("Nyquist Frequency: fr_c = %f\n", fr_c);
  printf("Time goes from 0 to %f\n\n",T);
  printf("Temporal series is %d data points long\n", nP);
  printf("Frequency series is %d data points long\n", nF);
 
  p = vector(1,nF);  fr = vector(1,nF);  
  period(time,x,nP, ofac,hifac, fr,p, nF, &nout, &jmax, &prob);
  //fasper(time,y,nP, ofac,hifac, fr,p, nF, &nout, &jmax, &prob);
  Norma_Total_a_J(p, 1, nF);
  Saving_to_File_float_float("pwsp", fr, p, nF, 0);
  printf("\nMin Frequency: (f_{1}, p) = (%f, %f),",fr[1],p[1]);
  printf("  Max Frequency: (f_{%d}, p) = (%f, %f)\n\n",nF,fr[nF],p[nF]);
  printf("Maximum Power Observed: (f_max, p) = (%f, %f)\n",fr[jmax],p[jmax]);
  printf("Significance of the maximum against the hypothesis of random noise:\n");
  printf("Probability = %f\n", prob);
  printf("First Power Spectrum has been estimated successfully!!\n");

  period(time,y,nP, ofac,hifac, fr,p, nF, &nout, &jmax, &prob);
  //fasper(time,y,nP, ofac,hifac, fr,p, nF, &nout, &jmax, &prob);
  Norma_Total_a_J(p, 1,nF);
  Saving_to_File_float_float("pwsp", fr, p, nF, 1);
  printf("\nMin Frequency: (f_{1}, p) = (%f, %f),",fr[1],p[1]);
  printf("  Max Frequency: (f_{%d}, p) = (%f, %f)\n\n",nF,fr[nF],p[nF]);
  printf("Maximum Power Observed: (f_max, p) = (%f, %f)\n",fr[jmax],p[jmax]);
  printf("Significance of the maximum against the hypothesis of random noise:\n");
  printf("Probability = %f\n", prob);
  printf("Second Power Spectrum has been estimated successfully!!\n");
  
  free_vector(fr,1,nF); free_vector(p,1,nF);
  free_vector(time, 1,No_of_Points);
  free_vector(x, 1,No_of_Points); 
  free_vector(y, 1,No_of_Points);
  printf("Power Spectrum has been estimated successfully!!\n");
   
  return (0);
}





