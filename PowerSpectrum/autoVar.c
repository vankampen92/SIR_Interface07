/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                          POWER SPECTRUM ANALYSIS                          */
/*	                                                                     */
/*                           David Alonso, 2004 (c)                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "pwsp.h"

void autocoVar(float *p, int n)
{
  /* The output is the autocorrelation function already stored in the
     same input vector p[]. 
  */
  int i, j, k;
  float mean, sum;
  float *a;

  sum =0.;
  for(i=1;i<=n;i++) sum += p[i];
  mean = sum/(float)n;
  
  a = vector(0,n-1);
  set_to_value_float(a, n, 0.);
  for(j=0; j<n/2; j++){ /* $j$ is the lag in the auto-covariance function */
    for(i=0; i<n-j; i++){
      a[j] += (p[1+i]-mean)*(p[1+i+j]-mean);   
    }
    a[j] = a[j]/(float)n;
  }
  
  for(k=0; k<n/2; k++){
    p[k+1] = a[k];
    p[n-k] = a[k];
  }

  free_vector(a,0,n-1);
}











