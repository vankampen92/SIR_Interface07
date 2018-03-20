#include "SIR.h"

void peak_to_Trough_Maxims_Minims(double *Time, int **Yp, 
				  int No, int Level,
				  double *Time_Max, double *Time_Min, 
				  double *Maxim, double *Minim, int *noVal)	   
{
  int i,j,k, Num;
  double t_M;
  int *Signum;
  float Delta, Delta_Max;

  Signum = ivector(0,No);
  set_to_value_double (Maxim, No, 0.); set_to_value_double (Minim, No, 0.);

  for(i=0; i<No; i++) 
    if(Yp[i][1] > Level) Signum[i] = 1;
    else                 Signum[i] =-1;
  
  Num = 0; i=1; Delta_Max = Yp[0][1] - Level; t_M = Time[0];
  while(i < No){
    if(Signum[i]*Signum[i-1] == 1){
      Delta = Yp[i][1] - Level;
      if(Signum[i] == 1){
        Delta_Max = MAX(Delta_Max, Delta); if(Delta == Delta_Max) t_M = Time[i]; 
      }
      else{
	Delta_Max = MIN(Delta_Max, Delta); if(Delta == Delta_Max) t_M = Time[i]; 
      }
      i++;
    }
    else{
      if(Signum[i-1] == 1){ Maxim[Num] = Delta_Max + Level; Time_Max[Num] = t_M; }
      else                { Minim[Num] = Delta_Max + Level; Time_Min[Num] = t_M; }
      Delta_Max = Yp[i][1] - Level; t_M = Time[i];
      k =i+1;
      while(k < No && Signum[k]*Signum[k-1] == 1){
	Delta = Yp[k][1] - Level;
	if(Signum[k] == 1){ 
	  Delta_Max = MAX(Delta_Max, Delta); if(Delta == Delta_Max) t_M = Time[k]; 
	}
	else{
	  Delta_Max = MIN(Delta_Max, Delta); if(Delta == Delta_Max) t_M = Time[k]; 
	}
	k++;
      }
      if(Signum[k-1] == 1){ Maxim[Num] = Delta_Max + Level; Time_Max[Num] = t_M; }
      else                { Minim[Num] = Delta_Max + Level; Time_Min[Num] = t_M; }
      Delta_Max = Yp[k][1] - Level; t_M = Time[k];
      i =k+1; 
      Num++;
    }
  }
  *noVal = Num;
  free_ivector(Signum,0, No);
}

float Relative_Amplification_Ratio(double *Maxim, int No, int I_0)
{
  int i;
  float ratio;
  float Valor;

  ratio = 0.;
  
  for(i=0; i < No; i++) ratio += ( Maxim[i] - (float)I_0 )/(float)I_0;
  
  ratio /= (float)No;
  return(ratio); 
}

float peak_to_Trough_Average_Ratio(double *Maxim, double *Minim, int No)
{
  int i;
  float Delta;
  
  Delta = 0.;
  for(i=0; i<No; i++){
    if(Maxim[i] != 0. && Minim[i] != 0.){
      Delta += (float)Maxim[i]/(float)Minim[i];
    }
  }
  Delta /= (float)No;
  return(Delta);
}

void saving_peak_to_Trough(double *Time_Max, double *Time_Min, 
			   double *Maxim, double *Minim, int No)
{
  int i;
  FILE *fp0, *fp1;

  fp0 = fopen("maxims.dat", "w"); fp1 = fopen("minims.dat", "w");
  
  for(i=0; i<No; i++){
    if(Maxim[i] != 0. && Minim[i] != 0.){
      printf("Maxim: %g\tTime: %g\n", Maxim[i], Time_Max[i]);
      printf("Minim: %g\tTime  %g\n", Minim[i], Time_Min[i]);
      fprintf(fp0, "%g\t%g\n", Time_Max[i], Maxim[i]);
      fprintf(fp1, "%g\t%g\n", Time_Min[i], Minim[i]);
    }
  }
  fclose(fp0); fclose(fp1);
}

float peak_to_Trough_Global_Range(double *Time, int **Yp, int n)
{
  int i;
  float ratio;
  float peak, trough;

  peak = (float)Yp[0][1]; trough = (float)Yp[0][1]; 
  i = 0;
  while(Yp[i][1] > 0 && i<n){
    trough = MIN(trough, Yp[i][1]);
    peak   = MAX(peak, Yp[i][1]);
    i++;
  }

  if(i == n){	
        if(Time[i-1] > 20.){ 
    	   ratio = peak/trough;
        }
	else{
	   ratio = 0.;
        }
  }
  else{
	if(Time[i] > 20.){
	   ratio = peak/trough;
        }
        else{	
    	   ratio = 0.;
	}
  }
  return (ratio);
}
  








