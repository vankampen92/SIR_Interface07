#include "../SIR.h" 
#include "../SIR_SEIR.h"

extern ParamSIR *Pi;
extern int TIME_DEPENDENT_PARAM; 
extern float timeFactor;

void derivaSIR (double x, double y[], double dydx[])
{ 
  float B;
  double Term, lnBeta;

  switch (TIME_DEPENDENT_PARAM) {
    case 0: /* Non seasonal forcing */
      B = Pi->b_0;
      break;	 
    case 1: /* Sinusoidal Seasonal Forcing */
      B = Pi->b_0 * (1. + Pi->b_1 * sin(2.* M_PI * x/Pi->Per));
      break;
    case 2: /* School term Seasonal Forcing (Keeling-Rohani-Grenfell, 2001) */
      Term = Time_of_the_Year(x/timeFactor);
      lnBeta = log(Pi->b_0) + Term * log(1. + Pi->b_1);   
      B = exp(lnBeta);
      break;
    case 3: /* School term Seasonal Forcing (Earn-Rohani-Bolker-Grenfell) */
      Term = Time_of_the_Year(x/timeFactor);
      B = Pi->b_0 * (1.+ Pi->b_1 * Term);
      break;
    default:
      printf("Invalid Seasonal Forcing (1, 2, 3)\n");
      printf(" TIME_DEPENDENT_PARAM = %d\n", TIME_DEPENDENT_PARAM);
      exit(0);
  }

  B = B/(y[0]+y[1]+y[2]); /* True mass action: Frequency dependene transition rate */

  dydx[0]  =  Pi->nu *(y[0]+y[1]+y[2]) + Pi->mu*y[2] - B * y[0]*y[1] - Pi->Imm *y[0] - Pi->d *y[0]; 
  dydx[1]  = +B * y[0]*y[1] + Pi->Imm *y[0] - Pi->g *y[1] - (Pi->d + Pi->a) *y[1]; 
  dydx[2]  =  Pi->g *y[1] - Pi->d *y[2] - Pi->mu *y[2]; 
}




