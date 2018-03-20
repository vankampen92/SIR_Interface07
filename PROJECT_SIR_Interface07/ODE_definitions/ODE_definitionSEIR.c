#include "../SIR.h" 
#include "../SIR_SEIR.h"

extern ParamSEIR *Pe;
extern int TIME_DEPENDENT_PARAM; 
extern float timeFactor;

void derivaSEIR (double x, double y[], double dydx[])
{ 
  float B;
  double Term, lnBeta;

  switch (TIME_DEPENDENT_PARAM) {
    case 0: /* Non seasonal forcing */
      B = Pe->b_0;
      break;	 
    case 1: /* Sinusoidal Seasonal Forcing */
      B = Pe->b_0 * (1. + Pe->b_1 * sin(2.* M_PI * x/Pe->Per));
      break;
    case 2: /* School term Seasonal Forcing (Keeling-Rohani-Grenfell, 2001) */
      Term = Time_of_the_Year(x/timeFactor);
      lnBeta = log(Pe->b_0) + Term * log(1. + Pe->b_1);   
      B = exp(lnBeta);
      break;
    case 3: /* School term Seasonal Forcing (Earn-Rohani-Bolker-Grenfell) */
      Term = Time_of_the_Year(x/timeFactor);
      B = Pe->b_0 * (1.+ Pe->b_1 * Term);
      break;
    default:
      printf("Invalid Seasonal Forcing (1, 2, 3)\n");
      printf(" TIME_DEPENDENT_PARAM = %d\n", TIME_DEPENDENT_PARAM);
      exit(0);
  }

  B = B/(y[0]+y[1]+y[2]+y[3]); /* True mass action: Frequency dependene transition rate */

  dydx[0]  =  Pe->nu * (y[0]+y[1]+y[2]+y[3]) + Pe->mu*y[3] - B * y[0]*y[2] - Pe->Imm *y[0] - Pe->d *y[0]; 
  dydx[1]  = +B * y[0]*y[2] + Pe->Imm*y[0] - Pe->si * y[1] - (Pe->d + Pe->a) *y[1]; 
  dydx[2]  =  Pe->si * y[1] - Pe->g *y[2] - (Pe->d + Pe->a) * y[2];
  dydx[3]  =  Pe->g *y[2] - Pe->d *y[3] - Pe->mu * y[3]; 
}




