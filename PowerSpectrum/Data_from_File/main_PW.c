/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	               Implementing spectral_Density(...)                    */ 
/*	            based on spctrm(...) form Numerical Recipes.             */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/*                      Data are read from a file introduced                 */
/*                              from the command line                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../pwsp.h"

        /* Global Shared parameters main Program <---> ArgumentControl() */
int m; /* Number of different frequencies in the Power Spectrum...       */
int k; /* Number of Segments */
int ovrlap;
char *File;
float time_Step;
int No_of_Data; 
int label;

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
  float Population[7]={8171000.,436000.,747500.,661000.,269400.,1105000.,494000.};
  FILE *fp;
  float *pe, *fr;
  float *X, *Xd, *pAcc;
  double Factor;
  int i, i_0, j;
  float time, Ave;
  char Fi[100] = "timeData.dat";

  /* Default parameter values */
  No_of_Data = 1000;
  File = Fi; time_Step = 7./365;  
  ovrlap = 1; m = 256; k = 3; label = 0;

  if(argc>1) 
    ArgumentControl(argc,argv); 

  if(ovrlap == 1)
    No_of_Points = (2*k + 1) * m;
  else 
    No_of_Points = 4*k*m;

  printf("Number of Data Points: %d\tNumber of used Points: %d\n", 
	 No_of_Data, No_of_Points);
  if(No_of_Data < No_of_Points){
    printf(" There are not enough data points!!\n  Try to choose a lower m.\n");
    printf(" Number of Frequences: %d\tNumber of Segments: %d\n", m, k);
    exit(0);
  }

  pAcc = vector(0, m); set_to_value_float(pAcc, m, 0.);  
  pe = vector(1,m+1);  fr = vector(0,m); 
  X  = vector(1,No_of_Data);
  Xd = vector(1, No_of_Points);

  fp = fopen(File, "r");
  for(i=1; i<=No_of_Data; i++){
    fscanf(fp, "%f\t%f\n", &time,&X[i]);
    //printf("Data(t = %f) = %f\n", time, X[i]);
  }
  fclose(fp);
  
  Ave = Average_float_Vector(X, No_of_Data);

  for(i=0; i<=m; i++) fr[i] = (float)i/2./(float)m/time_Step;
  
  Press_Key();
  printf("Transforming the data...\n");
  
  for(i_0=0; i_0<(No_of_Data-No_of_Points); i_0++){  
    for(i=1; i<=No_of_Points; i++){
      Xd[i] = (X[i_0+i]-Ave)/sqrt(Population[label]);
    }
    
    printf("... ...");  
    spectral_Density(Xd, pe, m, k, ovrlap);

    accummulating_ONE_PowerSpectrum(pe, pAcc, m); 
  }
  printf("\n");
  Factor = (double)(No_of_Data-No_of_Points);
  averaged_ONE_PowerSpectrum(Factor, pAcc, m, No_of_Data-No_of_Points, 0);

  Saving_to_File_float_float("pwsp", fr, pAcc, m, label, 0); 
  printf("Spectral Power per bin is saved\n\n");
  
  estimated_Spectral_Density(pAcc, fr, m);
  Saving_to_File_float_float("spden", fr, pAcc, m, label, 0);
  printf("Spectral Densitiy (power per unit frequence) is saved\n\n");

  free_vector(pe, 1,m+1); free_vector(pAcc, 0,m); 
  free_vector(fr, 0, m); 
  free_vector(Xd, 1, No_of_Points); free_vector(X, 1, No_of_Data); 
  printf("The power spectrum has been estimated successfully!!\n");
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}





