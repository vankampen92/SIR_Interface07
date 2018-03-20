#include "../SIR.h"
#include "../SIR_Analytic_General.h"

double Delta_2(ParamSet *P)
{
    double Discriminant;

    Discriminant = Gamma_2(P) - 2.*Omega_2(P);

    return Discriminant;
}

double Delta_4(ParamSet *P)
{
    double Discriminant;

    Discriminant = Gamma_2(P) - 4.*Omega_2(P);

    return Discriminant;
}

float Determinant(ParamSet *P, float b, float d)
{
  double De; 
  float D;
  ParamSet Pa;

  Pa.Beta = b; Pa.Delta = d;
  Pa.Alpha = P->Alpha; Pa.Mu = P->Mu; Pa.Imm = P->Imm; Pa.Gamma = P->Gamma;

  De = Omega_2(&Pa); D = (float)De;
  
  return (D);
}

float Discriminant_2(ParamSet *P, float b, float d)
{
  float Discriminant;
  ParamSet Pa;

  Pa.Beta = b; Pa.Delta = d;
  Pa.Alpha = P->Alpha; Pa.Mu = P->Mu; Pa.Imm = P->Imm; Pa.Gamma = P->Gamma;

  Discriminant = (float)Delta_2(&Pa);

  return Discriminant;
}

float Discriminant_4(ParamSet *P, float b, float d)
{
  float Discriminant;
  ParamSet Pa;

  Pa.Beta = b; Pa.Delta = d;
  Pa.Alpha = P->Alpha; Pa.Mu = P->Mu; Pa.Imm = P->Imm; Pa.Gamma = P->Gamma;

  Discriminant = (float)Delta_4(&Pa);

  return Discriminant;
}

float border_beta_delta(float (*Function)(ParamSet *, float, float), ParamSet *P, 
			float b0, float b1, float Tolerance, float d)
{
  /* 
     This function computes the \hat{beta} (aprox R_0) for any value of the relative infectious 
     phase, \hat{d}, corresponing to a value of zero for the 
     (*Function)(P, R_0, d) = Discriminant(P, R_0, d) = Tr^2 - 4.*Determinant = 0.
     where P is the parameter structure storing the values of the other model parameters.
  */ 
    float beta_R_0;
    float Zero; 
    int i;
    double F0, F1;
    

    F0 = (*Function)(P, b0,d); F1 = (*Function)(P, b1,d);

    if(F0*F1 > 0){
      printf("b = %g, a = %g, d = %g, mu = %g, g = %g, Imm = %g\n",
	     P->Beta, P->Alpha, P->Delta, P->Mu, P->Gamma, P->Imm);
      printf("Failure bracking the root\n");
      printf("Function(b=%f,d=%f) = %f\tFunction(b=%f,d=%f) = %f\n", b0,d,F0,b1,d,F1);  
      //exit(0);
      return (0.);
    }
   
    beta_R_0 = (double)zbrent_Param(Function, P, b0, b1, Tolerance, d);

    Zero = (*Function)(P, beta_R_0, d);
#if !defined SILENT
    printf("Checking for success in zero finding at %f level of tolerance...\n", Tolerance);
    printf("If a good root has been calculated, then: F(b, d) = 0., within certain tolerance\n");
    printf("F(%7.5f, %7.5f) = %f  +-  %f\n", beta_R_0, d, Zero, Tolerance); 
#endif

    return (beta_R_0);
}
