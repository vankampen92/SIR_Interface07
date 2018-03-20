#include "SIR.h"
#include "SIR_Analytic.h"

double Alpha(double a, double d, double g, double mu, double b, int I)
{
  double x;
  double a21,a12,a22,a11,b11,b12,b21,b22;
  
  a12 = a_12(a,d,g,mu,b);
  a21 = a_21(a,d,g,mu,b);
  a22 = a_22(a,d,g,mu,b);
  b11 = b_11(a,d,g,mu,b);
  b12 = b_12(a,d,g,mu,b);
  b22 = b_22(a,d,g,mu,b);

  if(I == 0){
    x = a21*a21*b11;                                    /* Eq. 2, for Susceptibles */
  }
  else if(I == 1){
    x = a12*a12*b22 - 2. *a22*a12*b12 + a22*a22*b11;    /* Eq. 1, for Infectives */
  }
  else{
    printf("Error in function Alpha()\n");
    return(0);
  }

  return (x);
}

double Beta(double a, double d, double g, double mu, double b, int I)
{
  double x;

  if(I == 0){                                            /* Eq. 2, for Susceptibles */
    x = b_22(a,d,g,mu,b);
  
  }
  else if(I == 1){
    x = b_11(a,d,g,mu,b);                                /* Eq. 1, for Susceptibles */
  }
  else{
    printf("Error in function Beta()\n");
    return(0);
  }

  return (x);
}

double Resonance_Frequency(double a, double d, double g, double mu, double b)
{
  double f;

  f = sqrt(Omega_2(a,d,g,mu,b)-0.5*Trace_2(a,d,g,mu,b))/2./M_PI;

  return (f);
}

double Peak_Frequency(double a, double d, double g, double mu, double b, int I)
{
  double f, A, x,y, wr;
  
  x = Alpha(a,d,g,mu,b,I);
  y = Beta(a,d,g,mu,b,I);   A = x/y;
  
  wr = 2.*M_PI * Resonance_Frequency(a,d,g,mu,b);

  f = 0.5*log(A*(A + 2.*wr*wr));
  
  f = -A + exp(f);
  
  f = sqrt(f)/2./M_PI;

  if(I == 0){
    printf("Alpha_S(%g,%g,%g,%g,%g) = %g\n",a,d,g,mu,b,x);
    printf("Beta_S(%g,%g,%g,%g,%g) = %g\n",a,d,g,mu,b,y);
  }
  else{
    printf("Alpha_I(%g,%g,%g,%g,%g) = %g\n",a,d,g,mu,b,x);
    printf("Beta_I(%g,%g,%g,%g,%g) = %g\n",a,d,g,mu,b,y);
  }
  return(f);
}

double Peak_Frequency_Approx(double a, double d, double g, double mu, double b, int I)
{
  double f, x,y, wr;
  
  y = Beta(a,d,g,mu,b,I);
  x = Alpha(a,d,g,mu,b,I);

  wr = 2.*M_PI * Resonance_Frequency(a,d,g,mu,b);

  f = wr*wr *(x - 0.5*y * wr*wr)/x;
  f = sqrt(f)/2./M_PI;

  return(f);
}

double Damping_Frequency(double a, double d, double g, double mu, double b)
{
  double f;

  f = sqrt(Omega_2(a,d,g,mu,b)-0.25*Trace_2(a,d,g,mu,b))/2./M_PI;

  return (f);
}
  
double Maximum_Power(double a, double d, double g, double mu, double b, int I)
{
  double p;
  double w_r;
  double x, GammaSqq;
  
  w_r = 2*M_PI*Resonance_Frequency(a,d,g,mu,b);
  GammaSqq = Trace_2(a,d,g,mu,b);
  x = (w_r*w_r - Omega_2(a,d,g,mu,b))*(w_r*w_r - Omega_2(a,d,g,mu,b)) + w_r*w_r*GammaSqq;
  p = (Alpha(a,d,g,mu,b,I) + Beta(a,d,g,mu,b,I) *w_r*w_r) /x;

  return(p);
}

double Delta_2(double a, double d, double g, double mu, double b)
{
    double Discriminant;

    Discriminant = Trace_2(a,d,g,mu,b) - 4.*Omega_2(a,d,g,mu,b);

    return Discriminant;
}

double Omega_2(double a, double d, double g, double mu, double b)
{ 
  double Det;
  
  Det = (d+mu)*(b-g-d-a);
  
  return (Det);
}

double Trace_2(double a, double d, double g, double mu, double b)
{ 
  double Trace; 

  Trace = -(d+mu)*(b+mu-a)/(g+d+mu);
  Trace *= Trace;

  return (Trace);
}

double Fi(double a, double d, double g, double mu, double b)
{
  /* Stable point for Infective */
  double I;

  I = (d+mu)/b*(b-g-d-a)/(g+d+mu);

  return (I);
}

double Psi(double a, double d, double g, double mu, double b)
{
  /* Stable Point for Susceptible */
  double S;

  S = (g+d+a)/b;

  return(S);
}

double a_12(double a, double d, double g, double mu, double b)
{
  double x;

  x = (d+mu)*(b-d-a-g)/(d+mu+g);

  return (x);
}

double a_21(double a, double d, double g, double mu, double b)
{
  double x;

  x = -(d+mu+g);

  return (x);
}
 
double a_22(double a, double d, double g, double mu, double b)
{
  double x;

  x = (d+mu)*(a-b-mu)/(d+mu+g);

  return (x);
}
 
double b_11(double a, double d, double g, double mu, double b)
{
  double x;
  
  x = b*Fi(a,d,g,mu,b)*Psi(a,d,g,mu,b) + (d+g+a)*Fi(a,d,g,mu,b);

  return (x);
}
 
double b_12(double a, double d, double g, double mu, double b)
{
  double x;

  x = -(d+a)*Fi(a,d,g,mu,b);

  return (x);
}

double b_21(double a, double d, double g, double mu, double b)
{
  double x;

  x = b_12(a,d,g,mu,b);

  return (x);
}
 
double b_22(double a, double d, double g, double mu, double b)
{
  double x;

  x = b*Fi(a,d,g,mu,b)*Psi(a,d,g,mu,b)+(d+mu)*(1.-Psi(a,d,g,mu,b))+(a-mu)*Fi(a,d,g,mu,b);

  return (x);
}








 









