#include "SIR.h"

int which_Village(double max_Probability, Community *Pop)
{
  double *a;
  int i,p;
  Community *P;

  a = dvector(0,No_of_Villages);
  P  = Pop;
  for(i=0; i<No_of_Villages; i++, P++)
    a[i] = P->ratePatch;
  
  p = Discrete_Sampling(a,No_of_Villages);
  //p = Discrete_Sampling_Rejection_Method(max_Probability, a, No_of_Villages);
  free_dvector(a, 0,No_of_Villages);
  
  return(p-1); 
}

struct point RandomPosition(int Length)
{
  struct point position;
  
  position.x = drand48()*Length;
  position.y = drand48()*Length;

  return position;
}

int Total_Number_of_Infective(int *SumPop, int Sp)
{
  int spNumber,k;

  spNumber = 0;
  for(k=1; k<=Sp; k++)
    if(SumPop[k] > 0) spNumber++;
  
  return spNumber;
}

void Fixed_Points(ParamSet *P, double *Psi, double *Fi)
{
  
  (*Fi) = (P->Delta+P->Mu)/P->Beta*(P->Beta-P->Gamma-P->Delta-P->Alpha)/(P->Gamma+P->Delta+P->Mu); /* Infective Stationary Level */
  (*Psi) = (P->Gamma+P->Delta+P->Alpha)/P->Beta; /* Susceptible Stationary Level */
  printf("Populations fractions at Equilibrium (Imm = 0)...\n");
  printf("Susceptible: %g\n", *Psi);
  printf("Infective:   %g\n", *Fi);
}

void Stability(ParamSet *P)
{
  double C, D, Det, Tra;
  double W_0, T_0;
  Det = (P->Delta+P->Mu)*(P->Beta-P->Gamma-P->Delta-P->Alpha);
  Tra = -(P->Delta+P->Mu)*(P->Beta+P->Mu-P->Alpha)/(P->Gamma+P->Delta+P->Mu);
  printf("Stability Analysis Isolated System (Imm = 0)...\n");
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
    printf("Characteristic Period of Damped oscillations, T_0 = %g\n", T_0);
  }  
  if(D > 0){
    W_0 = sqrt(Det -1/2.*Tra*Tra);
    W_0 /= (2.*M_PI);
    T_0 = 1./W_0;
    printf("Characteristic Frequency of Ressonant oscillations, w_0 = %g\n", W_0);
    printf("Characteristic Period of Ressonant oscillations, T_0 = %g\n", T_0);
  }  
}

float approximate_Endogenous_Period(ParamSet *P)
{
  float period;
  float D,R,L;

  D = 1./P->Gamma;  L = 1./P->Mu;  R = D*P->Beta;
  period = 2.*M_PI*sqrt(D*L/(R-1.));

  return(period);
}










