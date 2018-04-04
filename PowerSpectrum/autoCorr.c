/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	                                                                     */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "pwsp.h"

void autocorrel(float *p, int n)
{
  /* n must be multiple of 2 */
  /* The output is the autocorrelation function already stored in the
     same input vector p[]. 
  */
  int i, j, k, no;
  float mean, sum, var;
  float *a;

  assert(n%2 == 0);

  sum =0.;
  for(i=1;i<=n;i++) sum += p[i];
  mean = sum/(float)n;

  sum =0;
  for(i=1;i<=n;i++) sum += (p[i]-mean)*(p[i]-mean);
  var = sum/(float)n;
  
  a = vector(0,n/2-1);
  set_to_value_float(a, n/2, 0.);
  for(j=0; j<n/2; j++){
    no = 0;
    for(i=0; i<n-j; i++){
      no++;
      a[j] += (p[1+i]-mean)*(p[1+i+j]-mean);   
    }
    sum = (float)no * var;
    a[j] = a[j]/sum;
  }

  for(k=0; k<n/2; k++){
    if(k < (n/2 - 1)){
      p[k+1] = a[k];
      p[n-k] = a[k];
    }
    else{
      p[n/2] = a[k];
    } 
  }

  free_vector(a,0,n/2-1);
}











