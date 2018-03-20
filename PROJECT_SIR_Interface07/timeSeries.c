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

void Even_timeSeries(double *Time, int **Yp, int NoP, int nPop, 
		     float *time, float *x, float *y, int n, 
		     float t_0, float t_1)
{
  /* Building an even time series through interpolation */
  int i,k, no_k;
  float t, xx, dx;
  float *xa, *ya;
  
  xa = vector(1,n); ya = vector(1,n);
  k = 0;
  for(i=0; i<NoP; i++){
    if(Time[i] > t_0 && Time[i] < (t_1+0.5)){
      k++;
      x[k] = Yp[i][0]/(float)nPop;
      y[k] = Yp[i][1]/(float)nPop;
      time[k] = (float)Time[i];
    }
  }

  no_k = k;

  printf("Between %f and %f there are %d points: t[%d] = %f\n", t_0,0.5+t_1,k,k,time[k]);

  if(no_k < n){
    printf("Between %f and %f there are only %d points!\n", t_0,t_1+0.5,k);
    printf("Stochastic simulations should last longer: t_1 should be larger\n");
    exit(0);
  }
  if(time[no_k] < t_1){
    printf("Stochastic simulations should last shorter:t_1 should be smaller\n");
    exit(0);
  }
    
  for(i=1; i<=n; i++){ 
    t = t_0 + (float)(i-1)*(t_1 - t_0)/(float)n;
    polint(time, x, no_k, t, &xx, &dx);
    xa[i] = xx;
    polint(time, y, no_k, t, &xx, &dx);
    ya[i] = xx;
  }
  for(i=1; i<=n; i++){ 
    t = t_0 + (float)(i-1)*(t_1 - t_0)/(float)n;
    x[i] = xa[i];
    y[i] = ya[i];
    time[i] = t;
  } 
  free_vector(xa, 1,n); free_vector(ya, 1,n);
}

void pseudoEven_timeSeries(double *Time, int **Yp, int NoP, int nPop, 
			   float *time, float *x, float *y, int n, 
			   float t_0, float t_1)
{
  /* Building an pseudoeven time series through interpolation and selecting
     the Stationary Regim, i.e., between t_0 and t_1*/
  int i,k, no_k;
  float t, xx, dx;
  float *xa, *ya;
  
  xa = vector(1,n); ya = vector(1,n);
  k = 0;
  for(i=0; i<NoP; i++){
    if(Time[i] > t_0 && Time[i] < t_1){
      k++;
      x[k] = Yp[i][0]/(float)nPop;
      y[k] = Yp[i][1]/(float)nPop;
      time[k] = (float)Time[i];
    }
  }
  no_k = k;

  printf("Between %f and %f there are %d points: t[%d] = %f\n", t_0,t_1,k,k,time[k]);

  if(no_k < n){
    printf("Between %f and %f there are only %d points!\n", t_0,t_1,k);
    printf("Stochastic simulations should last longer: t_1 should be larger\n");
    exit(0);
  }  
    
  for(i=1; i<=n; i++){ 
    t = t_0 + (float)(i-1)*(time[no_k] - t_0)/(float)no_k;
    xa[i] = x[i];
    ya[i] = y[i];
  }

  free_vector(xa, 1,n); free_vector(ya, 1,n);
}

void timeSeries(double *Time, int **Yp, int NoP, int nPop, 
		float *time, float *x, float *y, int n, 
		float t_0, float t_1)
{
  /* 
     Building standarized time series of length $n$ between time $t_0$ and $t_1$
     at the stationary regim by using the community size and the relation 
     n/N = ave + x/sqrt(N); 
  */
  int i,k, no_k;
  float t, xx, dx;
  float *xa, *ya;
  float Ave_x, Ave_y;

  xa = vector(1,n); ya = vector(1,n);
  
  k = 0;
  for(i=0; i<NoP; i++){
    if(Time[i] > t_0 && Time[i] < t_1){
      k++;
      x[k] = (float)Yp[i][0];
      y[k] = (float)Yp[i][1];
      time[k] = (float)Time[i];
    }
  }
  no_k = k;

  Ave_x = Average_float_Vector(x, no_k); Ave_y = Average_float_Vector(y, no_k);

  /* Standarizing the time series to obtain standarized fluctuations... */
  k = 0;
  for(i=1; i<NoP; i++){
    if(Time[i] > t_0 && Time[i] < t_1){
	k++;
	x[k] = ((float)Yp[i][0] - Ave_x)/sqrt((double)nPop); 
	y[k] = ((float)Yp[i][1] - Ave_y)/sqrt((double)nPop); 
    }
  }
  
  printf("Between %f and %f there are %d points: t[%d] = %f\n", t_0,t_1, k, no_k, time[k]);

  if(no_k < n){
    printf("Between %f and %f there are only %d points!\n", t_0,t_1,k);
    printf("Stochastic simulations should last longer: t_1 should be larger\n");
    exit(0);
  }  
    
  for(i=1; i<=n; i++){ 
    t = t_0 + (float)(i-1)*(time[n] - t_0)/(float)n;
    xa[i] = x[i];
    ya[i] = y[i];
    time[i] = t;
  }

  free_vector(xa, 1,n); free_vector(ya, 1,n);
}













