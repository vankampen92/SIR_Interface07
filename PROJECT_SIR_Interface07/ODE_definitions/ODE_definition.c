#include "../SIR.h" 

extern double a;
extern double mu;
extern double b, b_0, b_1, Per;
extern double g;
extern double d;
extern double Imm;
extern int TIME_DEPENDENT_PARAM; 
extern float timeFactor;

void deriva (double x, double y[], double dydx[])
{ 
  float B;
  double Term, lnBeta;

  switch (TIME_DEPENDENT_PARAM) {
    case 0: /* Non seasonal forcing */
      B = b_0;
      break;	 
    case 1: /* Sinusoidal Seasonal Forcing */
      B = b_0 * (1. + b_1 * sin(2.* M_PI * x/Per));
      break;
    case 2: /* School term Seasonal Forcing (Keeling-Rohani-Grenfell, 2001) */
      Term = Time_of_the_Year(x/timeFactor);
      lnBeta = log(b_0) + Term * log(1. + b_1);   
      B = exp(lnBeta);
      break;
    case 3: /* School term Seasonal Forcing (Earn-Rohani-Bolker-Grenfell) */
      Term = Time_of_the_Year(x/timeFactor);
      B = b_0 * (1.+ b_1 * Term);
      break;
    default:
      printf("Invalid Seasonal Forcing (1, 2, 3)\n");
      printf(" TIME_DEPENDENT_PARAM = %d\n", TIME_DEPENDENT_PARAM);
      exit(0);
  }
  dydx[0]  = -B * y[0]*y[1] - Imm*y[0] + (d+mu) *(1. -y[0]-y[1]) + (d+a) *y[1]; 
  dydx[1]  = +B * y[0]*y[1] + Imm*y[0] - g *y[1] - (d+a) *y[1];  
}

