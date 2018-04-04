/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	               Implementing spectral_Density(...)                    */ 
/*	            based on spctrm(...) form Numerical Recipes.             */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/*                      Data are simulated and averaged over a               */ 
/* 			       number of simulations                         */
/*                          Sinusoidal Signal with noise                     */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../pwsp.h"

         /* Global Shared parameters main Program <---> ArgumentControl() */
int m;  /* Number of different frequencies in the Power Spectrum...       */
int k; /* Number of Segments */
int ovrlap;
char *File;
float time_Step;
int label;
double wr;
double noise;
double amplitud;
int REPLICAS;

int main(int argc, char **argv)
     /* 
	This main program calculates a spectral estimate of a data set by reading 
	it in as segments, windowning, Fourier transforming, and accumulating the 
	power spectrum. Data segments may or may not be overlapped at the decision of
	the user. 
     */	
{
  int No_of_Points;
  float fraction;  
  FILE *fp;
  float *pe, *fr;
  float *X, *Xd, *pAcc;
  double Factor;
  int i, i_0, j;
  float time, Ave;
  char Fi[100] = "timeData.dat";

  /* Default parameter values */
  REPLICAS = 100;
  File = Fi; time_Step = 7./365.;  
  ovrlap = 1; m = 256; k = 3; label = 0;
  amplitud = 2.5; noise = 0.5;  wr = 1.;

  if(argc>1) 
    ArgumentControl(argc,argv); 

  if(ovrlap == 1)
    No_of_Points = (2*k + 1) * m;
  else 
    No_of_Points = 4*k*m;

  printf("Number of used Points: %d\n", No_of_Points);

  pAcc = vector(0, m); set_to_value_float(pAcc, m, 0.);  
  pe = vector(1,m+1);  fr = vector(0,m); 
  X  = vector(1,No_of_Points);
  Xd = vector(1, No_of_Points);

  
    
  for(i_0=0; i_0<REPLICAS; i_0++){  
    /* Generating Pseudo Data */
    sinusoidal_Data(X, No_of_Points, time_Step, amplitud, wr, noise);
    
    /* End of Data Generations */
    
    Ave = Average_float_Vector(X, No_of_Points);
    for(i=1; i<=No_of_Points; i++){
      Xd[i] = (X[i]-Ave);
    }
    
    printf("... ...");  
    spectral_Density(Xd, pe, m, k, ovrlap);
    
    accummulating_ONE_PowerSpectrum(pe, pAcc, m); 
  }
  
  /* Time Domain */
  for(i=0; i<=m; i++) fr[i] = (float)i * time_Step;
  Saving_to_File_float_float("serie", fr, X, m, label, 0);
  /* Frequence Domain */  
  for(i=0; i<=m; i++) fr[i] = (float)i/2./(float)m/time_Step;
  
  printf("\n");
  Factor = (double)(REPLICAS);
  averaged_ONE_PowerSpectrum(Factor, pAcc, m, No_of_Points, 0);

  Saving_to_File_float_float("pwsp", fr, pAcc, m, label, 0); 
  printf("Spectral Power per bin is saved\n\n");
  
  estimated_Spectral_Density(pAcc, fr, m);
  Saving_to_File_float_float("spden", fr, pAcc, m, label, 0);
  printf("Spectral Densitiy (power per unit frequence) is saved\n\n");

  free_vector(pe, 1,m+1); free_vector(pAcc, 0,m); 
  free_vector(fr, 0, m); 
  free_vector(Xd, 1, No_of_Points); free_vector(X, 1, No_of_Points); 
  printf("The power spectrum has been estimated successfully!!\n");
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}





