#include "../SIR.h"
#include "../SIR_Analytic_General.h"

float resilience(ParamSet *P)
{ 
  float resil;
  
  resil = -Real_Eigen_Value(P);
  
  if(resil < 0.) resil = 0.; /* The system is not resilient */
  
  return (resil);
}

float Real_Eigen_Value(ParamSet *P)
{
  /* This function return the real part of the eigen values,
     if they are complex numbers or the maximum value of them 
     if both are real; */

  int stability_Condition, damping_Condition;
  float lambda_1, lambda_2;
  double Re_lambda, Im_lambda;
  double Det, Tra; 

  Det = Omega_2(P);
  Tra = Trace(P);  
  stability_Condition = Condition_Stability(P);

  if(stability_Condition == 1){
    damping_Condition = Condition(P,0.25);
    if(damping_Condition == 1){
      Re_lambda = 0.5 *Tra;
      Im_lambda = 0.5 *sqrt(4.*Det - Tra*Tra);
      lambda_1 = Re_lambda;
    }
    else{
      lambda_1 = 0.5 *(Tra + sqrt(Tra*Tra - 4.*Det));
      lambda_2 = 0.5 *(Tra - sqrt(Tra*Tra - 4.*Det));
      lambda_1 = MAX(lambda_1,lambda_2);
    }
  }
  else{
    lambda_1 = 0.5 *(Tra + sqrt(Tra*Tra - 4.*Det));
    lambda_2 = 0.5 *(Tra - sqrt(Tra*Tra - 4.*Det));
    lambda_1 = MAX(lambda_1,lambda_2);
  }
  return(lambda_1);
}

float Imaginary_Eigen_Value(ParamSet *P)
{
  /* This function return the real part of the eigen values,
     if they are complex numbers or the maximum value of them 
     if both are real; */

  int stability_Condition, damping_Condition;
  float lambda_1, lambda_2;
  double Re_lambda, Im_lambda;
  double Det, Tra; 

  Det = Omega_2(P);
  Tra = Trace(P);  
  stability_Condition = Condition_Stability(P);

  if(stability_Condition == 1){
    damping_Condition = Condition(P,0.25);
    if(damping_Condition == 1){
      Re_lambda = 0.5 *Tra;
      Im_lambda = 0.5 *sqrt(4.*Det - Tra*Tra);
      lambda_1 = Im_lambda;
    }
    else{
      lambda_1 = 0.;
    }
  }
  else{
    lambda_1 = 0.;
  }
  return(lambda_1);
}
