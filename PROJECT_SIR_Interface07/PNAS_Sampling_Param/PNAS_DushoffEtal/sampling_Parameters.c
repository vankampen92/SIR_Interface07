#include "../../SIR.h"

extern double b_0;

void sample_Parameters(ParamSet *P, 
		       float *range_R, float *range_L, float *range_D)
{
  float R,L,D;
  float range;
  /* Sampling R_0 */
  range = range_R[1] - range_R[0];
  R = range_R[0] + drand48()*range;
  /* Sampling L */
  range = range_L[1] - range_L[0];
  L = range_L[0] + drand48()*range;
  /* Sampling D */
  range = range_D[1] - range_D[0];
  D = range_D[0] + drand48()*range;
   
  P->Mu = 1./L;
  P->Gamma = 1./D;
  b_0 = P->Beta = R * P->Gamma;
  printf("Transmission Rate: %g\tAverage Immunity Period: %f y.\tAverage Infectious Period: %f d.\n", P->Beta, L, 365.*D);  
}
  








