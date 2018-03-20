#include "SIR.h"
#include "SIR_Analytic.h"

float Delta_2_DimLess_4(float b, float d)
{
    float Discriminant;
    double a,g,mu;

    mu = 0.;  a = 0.;  g = 1.;
 
    Discriminant = Trace_2(a,(double)d,g,mu,(double)b) - 4.*Omega_2(a,(double)d,g,mu,(double)b);

    return Discriminant;
}

float Delta_2_DimLess_2(float b, float d)
{
    float Discriminant;
    double a,g,mu;

    mu = 0.;  a = 0.;  g = 1.;
 
    Discriminant = Trace_2(a,(double)d,g,mu,(double)b) - 2.*Omega_2(a,(double)d,g,mu,(double)b);

    return Discriminant;
}

float Discriminant_Zero_2(float b, float d)
{
  float Discriminant;

  Discriminant = d*b*b/(1.+d)/(1.+d) - 2.*(b-d-1);

  return Discriminant;
}

float Discriminant_Zero_4(float b, float d)
{
  float Discriminant;

  Discriminant = d*b*b/(1.+d)/(1.+d) - 4.*(b-d-1);

  return Discriminant;
}

float beta_R_0_delta(float (*Function)(float, float), 
		     float b0, float b1, float Tolerance, float d)
{
  /* 
     This function computes the beta (or R_0) for any value of the relative infectious 
     phase length, d corresponing to a value of zero for the 
     (*Function)(R_0, d) = Discriminant(R_0, d) = Tr^2 - 4. Determinant = 0. 
  */ 
    float beta_R_0;
    float Zero; 
    int i;
    double F0, F1;
    float zbrent_1(float (*)(float, float), float, float, float, float);

    F0 = (*Function)(b0,d); F1 = (*Function)(b1,d);
  
    if(F0*F1 > 0){
	printf("Failure bracking the root\n");
	exit(0);
    }
   
    beta_R_0 = (double)zbrent_1(Function, b0, b1, Tolerance, d);

    Zero = (*Function)(beta_R_0, d);
#if !defined SILENT
    printf("Checking for success in zero finding at %f level of tolerance...\n", Tolerance);
    printf("If beta_R_0 is a good root then: F(b, d) = 0., within certain tolerance\n");
    printf("F(%7.5f, %7.5f) = %f  +-  %f\n", beta_R_0, d, Zero, Tolerance); 
#endif

    return (beta_R_0);
}










