/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "SIR.h"

void timeSeries_vanKampen(double Ave_x, double Ave_y,
			  double *Time, int **Yp, int NoP, int nPop, 
			  float *time, float *x, float *y, int n, 
			  float t_0, float t_1)
{
  /* 
     Building standarized time series of length $n$ between time $t_0$ and $t_1$
     at the stationary regim by using the community size and the relation 
     $n/N = phi + x/sqrt(N)$, where $phi$ is the stationary value at the 
     equilibrium of the corresponding deterministic system; 
  */
  int i,k, no_k;
  float t, xx, dx;
  float *xa, *ya;

  k = 0;
  for(i=0; i<NoP; i++){
    if(Time[i] > t_0 && Time[i] < t_1){
      k++;
      time[k] = (float)Time[i];
      x[k] = ((float)Yp[i][0] - (float)nPop*(float)Ave_x)/sqrt((double)nPop); 
      y[k] = ((float)Yp[i][1] - (float)nPop*(float)Ave_y)/sqrt((double)nPop); 
    }
  }
  no_k = k;
  
  printf("Between %f and %f there are %d points: t[%d] = %f\n", t_0,t_1, k, no_k, time[k]);
  if(no_k < n){
    printf("Between %f and %f there are only %d points!\n", t_0,t_1,k);
    printf("Stochastic simulations should last longer, i.e, t_1 should be larger\n");
    exit(0);
  }  
    
  for(i=1; i<=n; i++){ 
    t = t_0 + (float)(i-1)*(time[n] - t_0)/(float)n;
    time[i] = t;
  }
}













