#include "../SIR.h"
#include "../SIR_Analytic_General.h"

/* Stable point in the ODE system (Fi, Psi), Fi being the fraction 
   of the population described by the first equation (usually,
   fraction of susceptibles), and Psi the fraction of the population
   described by the second equation (usually, fraction of infective) 
   of our ODE system. 
*/
double Fi(ParamSet *P)
{
  /* Stable Point for Susceptible */
  double S, signe;
  double a, b, c;
 
  if( P->Beta >= (P->Delta + P->Gamma + P->Alpha) )
    signe = -1;
  else
    signe = +1;
  
  a = 1.; 
  b = -(1.+ (P->Delta+P->Gamma+P->Alpha)/P->Beta + (P->Delta+P->Gamma+P->Mu)/(P->Delta+P->Mu)/P->Beta *P->Imm);
  c = (P->Delta+P->Alpha+P->Gamma)/P->Beta;

  S = 0.5*(-b + signe * sqrt(b*b - 4.*a*c))/a;

  return(S);
}
double Psi(ParamSet *P)
{
  /* Stable point for Infective */
  double I;

  I = (P->Delta+P->Mu)*(1.-Fi(P))/(P->Delta+P->Gamma+P->Mu);

  return (I);
}
void Fixed_Points_General(ParamSet *P, double *S, double *I)
{
  (*S) = Fi(P);
  (*I) = Psi(P);
  printf("Populations fractions at Equilibrium (Imm = %g)...\n", P->Imm);
  printf("Susceptible: %g\n", *S);
  printf("Infective:   %g\n", *I);
}
/* End of Fix Point expressions */

/* Entries of the stability (Jacobian) matrix A = (a_{ij}) 
   evalutated at the fix point (Fi, Psi) */
double a_11(ParamSet *P)
{
    /* Entry of a Stability matrix evaluated at the stable point of interest */
    double x;
    
    x = -P->Beta*Psi(P) -P->Imm -P->Delta -P->Mu;
    
    return (x);
}

double a_12(ParamSet *P)
{
    /* Entry of the  Stability matrix evaluated at the stable point of interest */
    double x;

    x = -P->Beta*Fi(P) +P->Alpha -P->Mu;

    return (x);
}

double a_21(ParamSet *P)
{
    /* Entry of a Stability matrix evaluated at the stable point of interest */
    double x;
    
    x = P->Beta*Psi(P) + P->Imm;
    
    return (x);
}
 
double a_22(ParamSet *P)
{
    /* Entry of a Stability matrix evaluated at the stable point of interest */
    double x;
    
    x = +P->Beta*Fi(P) - (P->Gamma+P->Delta+P->Alpha);
    
    return (x);
}

/* Entries of the noise covariance matrix B = (b_{ij}): 
   b_11 = alpha_{2,0}(Fi,Psi); b_22 = beta_{2,0}(Fi,Psi), b_12=b_21= gamma(Fi,Psi) 
   in the van Kampen (1992, chap 10) Large N expansion.
*/
double b_11(ParamSet *P)
{    
    /* b_11 = alpha_{2,0}(Fi,Psi)  */
    double x;
  
    x = P->Beta*Fi(P)*Psi(P) + P->Imm*Fi(P) + (P->Delta+P->Mu)*(1.-Fi(P)) + (P->Alpha-P->Mu)*Psi(P);
    
    return (x);
}
 
double b_12(ParamSet *P)
{
    /*  b_12 = gamma(Fi, Psi) */
    double x;
    
    x = -P->Beta*Fi(P)*Psi(P) -(P->Delta+P->Alpha)*Psi(P) -P->Imm*Fi(P);
    
    return (x);
}

double b_21(ParamSet *P)
{
    double x;
    
    x = b_12(P);
    
    return (x);
}

double b_22(ParamSet *P)
{
  /* b_22 = beta_{2,0}(Fi,Psi) */
  double x;
  
  x = P->Beta*Fi(P)*Psi(P) + (P->Delta+P->Alpha+P->Gamma)*Psi(P) + P->Imm*Fi(P);
  
  return (x);
}

double Omega_2(ParamSet *P)
{ 
  double Determinant;
  
  Determinant = a_11(P)*a_22(P) - a_12(P)*a_21(P);
  
  return (Determinant);
}

double Gamma_2(ParamSet *P)
{ 
  /* Gamma_2 = Trace_2*Trace_2 */
  double Trace_2; 

  Trace_2 = (a_11(P) + a_22(P)) * (a_11(P) + a_22(P));

  return (Trace_2);
}

double Trace(ParamSet *P)
{ 
  double Tr; 

  Tr = a_11(P) + a_22(P);

  return (Tr);
}

int Condition_Stability(ParamSet *P)
{
  int bool;
  
  bool = 0;
  if(Omega_2(P)>0. && Trace(P)<0.){
      bool = 1;
  }

  return (bool);
}

int Condition(ParamSet *P, double factor)
{
  int bool;
  double Discriminant;
  
  bool = 0;
  if(Omega_2(P)>0. && Trace(P)<0.){
    Discriminant = Omega_2(P) - factor*Gamma_2(P);
    if(Discriminant > 0.) 
      bool = 1;
  }

  return (bool);
}

void eigen_Values(ParamSet *P)
{
  double lambda_1, lambda_2;
  double Re_lambda, Im_lambda;
  double Det, Tra;
  int stability_Condition, damping_Condition;

  Det = Omega_2(P);
  Tra = Trace(P);
  
  stability_Condition = Condition_Stability(P);

  if(stability_Condition == 1){
    damping_Condition = Condition(P,0.25);
    if(damping_Condition == 1){
      Re_lambda = 0.5 *Tra;
      Im_lambda = 0.5 *sqrt(4.*Det - Tra*Tra);
      printf("lambda_1 = %g + i %g\n", Re_lambda, Im_lambda);
      printf("lambda_2 = %g - i %g\n", Re_lambda, Im_lambda);
    }
    else{
      lambda_1 = 0.5 *(Tra + sqrt(Tra*Tra - 4.*Det));
      lambda_2 = 0.5 *(Tra - sqrt(Tra*Tra - 4.*Det));
      printf("lambda_1 = %g\n", lambda_1);
      printf("lambda_2 = %g\n", lambda_2);

    }
  }
  else{
    lambda_1 = 0.5 *(Tra + sqrt(Tra*Tra - 4.*Det));
    lambda_2 = 0.5 *(Tra - sqrt(Tra*Tra - 4.*Det));
    printf("lambda_1 = %g\n", lambda_1);
    printf("lambda_2 = %g\n", lambda_2);
  } 
}

void Stability_General(ParamSet *P)
{
  double C, D, Det, Tra;
  double W_0, T_0;
  Det = Omega_2(P);
  Tra = Trace(P);
  printf("Stability Analysis Non Isolated SIR System (Imm = %g)...\n", P->Imm);
  printf("Det(A) = %g\tTra(A) = %g\n", Det, Tra);
  printf("Fix Point is Stable when Det(A) > 0 and Tra(A) < 0\n");
  printf("Stable Point shows damped oscillations if C = Det(A)-1/4*Tra^2(A) > 0\n");
  printf("Stable Point shows internally-induced resonant cycles if D = Det(A)-1/2*Tra^2(A)>0\n");
  C = Det -1/4.*Tra*Tra; D = Det -1/2.*Tra*Tra; 
  printf("C = %g\t",C); printf("D = %g\n",D);
  if(C > 0){
    W_0 = sqrt(Det -1/4.*Tra*Tra);
    W_0 /= (2.*M_PI);
    T_0 = 1./W_0;
    printf("Characteristic Frequency of Damped oscillations, w_0 = %g\n", W_0);
    printf("Characteristic Period of Damped oscillations, T_0 = %g = %g y.\n", T_0, T_0/365.);
  }  
  if(D > 0){
    W_0 = sqrt(Det -1/2.*Tra*Tra);
    W_0 /= (2.*M_PI);
    T_0 = 1./W_0;
    printf("Characteristic Frequency of Ressonant oscillations, w_0 = %g\n", W_0);
    printf("Characteristic Period of Ressonant oscillations, T_0 = %g = %g y.\n\n", 
	   T_0, T_0/365.);
  }  
}











