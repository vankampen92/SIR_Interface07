#include "../SIR.h"
#include "../SIR_Analytic_General.h"

double power_Spectrum(ParamSet *P, int I, double w)
{
    double p;
    /* This power spectrum is two-sided. Therefore you need a factor of two 
       Furthermore, when w < 0, the value will be considered to be zero */
    
    if(w >= 0)
      p = 2.*(Alpha(P,I) + Beta(P,I)*w*w)/((w*w - Omega_2(P))*(w*w - Omega_2(P)) + w*w* Gamma_2(P)); 
    else
      p = 0.;

    return p;
}

double Alpha(ParamSet *P, int I)
{
  /* I = 0, Fist equation; I = 1, Second Equation in your ODE system */
  double x;
  double a11,a21,a12,a22,b11,b12,b21,b22;
  
  a11 = a_11(P); a12 = a_12(P);
  a21 = a_21(P); a22 = a_22(P);
  b11 = b_11(P);
  b12 = b_12(P);
  b22 = b_22(P);

  if(I == 0){
    x = a12*a12*b22 - 2. *a22*a12*b12 + a22*a22*b11;
  }
  else if(I == 1){
    x = a21*a21*b11 - 2. *a11*a21*b12 + a11*a11*b22;
  }
  else{
    printf("Error in function Alpha()\n");
    exit(0);
  }

  return (x);
}

double Beta(ParamSet *P, int I)
{
  /* I = 0, Fist equation; I = 1, Second Equation in your ODE system */
  double x;

  if(I == 0){
    x = b_11(P);
  
  }
  else if(I == 1){
    x = b_22(P);
  }
  else{
    printf("Error in function Beta()\n");
    exit(0);
  }

  return (x);
}

double Resonance_Frequency(ParamSet *P)
{
  /* McKane-Newman Approx. */
  double f;

  f = sqrt(Omega_2(P)-0.5*Gamma_2(P))/2./M_PI;

  return (f);
}

double Resonance_Frequency_Peak_Approx(ParamSet *P, int I)
{
  /* Approximated Result */
  double f, x, e, wr;
  
  wr = 2.*M_PI*Resonance_Frequency(P); 
  x = Alpha(P,I)/Beta(P,I);
  e = Omega_2(P)/x; 

  f = x*(-1. + sqrt(1+e*e) + 1./x * 1./(sqrt(1.+e*e)) * wr*wr);
  f = sqrt(f)/2./M_PI;

  return(f);
}

double Damping_Frequency(ParamSet *P)
{
  double f;

  f = sqrt(Omega_2(P)-0.25*Gamma_2(P))/2./M_PI;

  return (f);
}
  
double Maximum_Power(ParamSet *P, int I)
{
  double p;
  double w_r;
  double x, GammaSqq;
  
  w_r = 2*M_PI*Resonance_Frequency_Peak(P, I);
 
  p = power_Spectrum(P, I, w_r);

  return(p);
}

double Resonance_Frequency_Peak(ParamSet *P, int I)
{
  /* Exact Result */
  double f, r, epsi;
  double a, b;

  r = Alpha(P,I)/Beta(P,I); assert(r > 0.);
  a = Omega_2(P); b = Gamma_2(P);  epsi = a/r;
   
  f = r*(-1. + sqrt(1.+ epsi*epsi + 2.*epsi - b/r));

  assert(f > 0.);
  f = sqrt(f)/2./M_PI;
  return (f);
}

int Exact_Condition_Peak(ParamSet *P, int I)
{
  int bool;
  double Discriminant;
  double r, epsi;
  double a, b;

  bool = 0; 
  if(Omega_2(P)>0. && Trace(P)<0.){
    r = Alpha(P,I)/Beta(P,I); assert(r > 0.);
    a = Omega_2(P); b = Gamma_2(P);  epsi = a/r;
   
    Discriminant = epsi*epsi + 2.*epsi - b/r;
    if(Discriminant > 0.) bool = 1;
  }

  return (bool);
}

/* double overall_Analytic(ParamSet *P, int I) */
/* { */
/*   double amplification, omega_2, gamma, gamma_2, A; */

/*   omega_2 = Omega_2(P); */
/*   gamma_2 = Gamma_2(P); */
/*   A = sqrt(omega_2*gamma_2*(4.*omega_2 - gamma_2)); */

/*   amplification = 2.* M_PI *(Alpha(P,I) + Beta(P,I)*omega_2)/A; */

/*   return(amplification); */
/* } */

double overall_Analytic(ParamSet *P, int I)
{
  double amplification, omega_2, gamma_2, A;

  omega_2 = Omega_2(P);
  gamma_2 = Gamma_2(P);
  A = sqrt(gamma_2)*omega_2;

  amplification = M_PI *(Alpha(P,I) + Beta(P,I)*omega_2)/A;

  return(amplification);
}

double coherence_Analytic(ParamSet *P, int I, double f_semi, double nu_p)
{
  /* nu_p is the dominant frequency: frequency at which the PSD peaks 
     in T^{-1}. Other frequencies here are all angular frequencies */
  double coherence, O_2, G_2, G, a, w1, w2, w_1p, w_1n, w_2p, w_2n;
  
  O_2 = Omega_2(P);
  G_2 = Gamma_2(P);
  G   = sqrt(G_2);
  a = sqrt(4.*O_2-G_2)/2.;
  
  w2 = 2.*M_PI*(nu_p+0.5*f_semi);
  w1 = 2.*M_PI*(nu_p-0.5*f_semi);  

  if(w1 < 0.){
    printf("Frequency band is too broad\n");
    printf("Coherence is calculated for frequency band... (%g, %g)\n",
	   nu_p-0.5*f_semi, nu_p+0.5*f_semi);
    exit(0);
  }
 
  w_1p = 2.*(w1+a)/G;     w_1n = 2.*(w1-a)/G;
  w_2p = 2.*(w2+a)/G;     w_2n = 2.*(w2-a)/G;
  
  coherence = (Alpha(P,I) - Beta(P,I)*O_2)/(4.*a*O_2) * log((w_2p*w_2p + 1.)*(w_1n*w_1n + 1.)/(w_2n*w_2n + 1.)/(w_1p*w_1p + 1.)) + (Alpha(P,I) + Beta(P,I)*O_2)/(G*O_2)*(atan(w_2p)+atan(w_2n)-atan(w_1p)-atan(w_1n));
  
  return(coherence);
}












