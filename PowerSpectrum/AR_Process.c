/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        AUTOREGRESSIVE MODELS OF ORDER p                   */
/*	                                                                     */
/*                             David Alonso, 2004 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "pwsp.h"

void AR_p(float *X, int n, int p, float sigma_2, float *alpha)
{
  /* AR process of simga^2 = 1, and coeficients all equal to 0.5 */ 
  int i, k;
  float AR;
  float y;

  /* Initial Conditions */
  for(i=1; i<=p; i++) {
    X[i] = sqrt(sigma_2) * gasdev_01();
  }
  
  /* AR(p) process */
  for(i=3; i<=n; i++){
    AR = 0.;
    for(k=1; k<=p; k++) AR += alpha[k] * X[i-k];
    X[i] = AR + sqrt(sigma_2) * gasdev_01();
  }
}

void AR_p_PWSP_Analytical(float *qu, int n, int p, 
			  float sigma, float *a, 
			  float fraction)
{
  int i;
  float f;
  
  if(p == 2){
    for(i=0; i<n; i++){
      f = (float)i * fraction/(float)n; 
      qu[i] = sigma/(1.+ a[1]*a[1] + a[2]*a[2] + 2.*a[2] + 2.*(a[1]*a[2]-a[1])*cos(2.*M_PI*f) -
		     4.*a[2]*cos(2.*M_PI*f)*cos(2.*M_PI*f))/(2.*M_PI); 
    }
  }
  else{
    printf("AR(%d) is not of order p = 2\n", p);
    printf("The analytical Spectrum cannot be calculated\n");
  }
}
