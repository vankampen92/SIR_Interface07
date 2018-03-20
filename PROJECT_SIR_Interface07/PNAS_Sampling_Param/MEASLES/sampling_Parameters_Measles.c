#include "../../SIR.h"

extern double b_0;

void sample_Parameters(ParamSet *P, 
		       float *range_L, float *range_D, float *range_R)
{
  /* Dimensionless Parameters Ranges:
     l_0 = 10.;                   l_1 = 20.;                 Range in Betas (transmission rates) 
     d_0 = 6.0e-4;                d_1= 2.0e-3.;              Range in Deltas (turnover rates) 
     r_0 = 1.e-5/time_Factor;     r_1 = 1.e-4/time_Factor;   Range in Immigration rates 
  */
  float R,L,D;
  float range;
  /* Sampling Betas */
  range = range_L[1] - range_L[0];
  L = range_L[0] + drand48()*range;
  /* Sampling Deltas */
  range = range_D[1] - range_D[0];
  D = range_D[0] + drand48()*range;
  /* Sampling Immigration rates */
  range = range_R[1] - range_R[0];
  R = range_R[0] + drand48()*range;

  P->Beta = (double)L; b_0 = P->Beta;
  P->Delta = (double)D;
  P->Imm = (double)R;
  
  printf("Transmission Rate: %g\t: Turnover rate (Death-Birth): %g\tImmigration %g\n",
	 P->Beta, P->Delta, P->Imm);
  printf("Rates in dimensionless paremeters\n");
}
  








