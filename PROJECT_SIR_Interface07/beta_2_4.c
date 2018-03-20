#include "SIR.h"
#include "SIR_Analytic.h"

float beta_2(float d, float signe)
{
  /* Solutions of Tr^2(M(d,beta) - 2.*Det(M(d,beta)) = 0., where M(d,beta) is the 
     stability matrix evaluated at the fix point */
  float beta_2;

  beta_2 = (1.+d)*(1.+d)/d * (1. + signe * sqrt((1.-d)/(1.+d)));

  return beta_2;
}

float beta_4(float d, float signe)
{
  /* Solutions of Tr^2(M(d,beta) - 4.*Det(M(d,beta)) = 0., where M(d,beta) is the 
     stability matrix evaluated at the fix point */
  float beta_4;

  beta_4 = 2.*(1.+d)*(1.+d)/d * (1. + signe * sqrt(1./(1+d)));

  return beta_4;
}
 









