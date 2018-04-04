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

void sinusoidal_Data(float *X, int No, 
		     float Step, double amplitud, double wr, double noise)
{
  int i;
  double time;

  for(i=1; i<=No; i++){
    time = (double)i*(double)Step;
    X[i] = amplitud * (sin(2.*M_PI*wr*time) + drand48()*noise);
    //printf("X[%d] = %f ", i, X[i]);
  }

}

