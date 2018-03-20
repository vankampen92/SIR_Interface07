#include "../../SIR.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern double a;
extern double mu;
extern double b;
extern double g;
extern double d;
extern double Imm;
extern int POPULATION;
extern int No_of_Points;
extern double factor;
extern float r;

void default_values()
{
  /* Initial settings and default values (Measles Values as a default) * * * * * * * *  */

  a = 0.;
  mu = 0.;
  b = 1.175;
  g = 0.077;
  d = 5.5e-5;
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  POPULATION = 100000;
  No_of_Points =  1000;
  factor = 5.; /* The analytic spectrum is calculated form f=0 to f=factor*f_M, where
		  f_M is the frequency at which the power spectrum peaks */
  r = 0.05;
  /* END (Initial settings and default values) * * * * * * * * * * * * * * * * */
}




