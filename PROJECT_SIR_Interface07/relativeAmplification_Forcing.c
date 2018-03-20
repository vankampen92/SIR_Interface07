#include "SIR.h"

extern double mu;
extern double b;
extern double b_0, b_1;
extern double g;
extern double Imm;
extern double a;
extern double d;

float Relative_Amplification_Ratio_Forcing(double *Time_Max, double *Maxim, int No,
					   int N, double f_0, double f_1,
					   ParamSet *Par, int No_of_Points)
{
  int i;
  float x_i, x_s, ratio;
  float *P;
  float Valor;

  P = vector(0,1);
  P[0] = f_0; P[1] = f_1;
  b_0 = b = Par->Beta;  g = Par->Gamma; mu = Par->Mu; a = Par->Alpha;
  Imm = Par->Imm; d = Par->Delta;
  x_i = 0.; /* Initial time */
  ratio = 0.;
  for(i=0; i < No; i++){
    x_s = Time_Max[i];
    integration(P, 2, x_i, x_s, No_of_Points, deriva);  
    x_i = x_s;
    Valor = P[1] * (float)N;
    ratio += fabs(Maxim[i]-Valor) / Valor;
    //printf("- Time: %g\tStochastic Value: %g\tDeterministic Value: %f\n", 
    //       Time_Max[i], Maxim[i], Valor);
  }
  ratio /= (float)No;
  
  printf("Successful numerical integration!\n");
  free_vector(P, 0,1);
  return(ratio);
}




  








