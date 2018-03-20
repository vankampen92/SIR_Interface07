#include "../../../SIR.h"
#include "../../../SIR_Analytic_General.h"
#define ITMAX 100
#define EPS 3.0e-8

float endogen_Period(ParamSet *P, float b, float d, float Period, int SI)
{
  float Value;
  double f_peak;

  P->Beta = b; P->Delta = d;
  
  f_peak = (float)Resonance_Frequency_Peak(P, SI) * 365 / 13;
  Value = 1/f_peak - Period; 

  return Value;
}

float border_Period(float (*Function)(ParamSet *, float, float, float, int), 
		    ParamSet *P, float b0, float b1, float Pe, int SI,
		    float Tolerance, float d)
{
  /* 
     This function computes the relative Beta (b/g) for any value of the relative infectious 
     phase, d/g, corresponing to a value of zero of the function:
     
     float (*Function)(P, Beta, Delta, Pe, SI)
     
     where P is the parameter structure storing the values of the other model parameters.
  */ 
    float beta_R_0;
    float Zero; 
    int i;
   
    beta_R_0 = (double)zbrent_Endogen(Function, P, b0, b1, Pe, SI, 
				      Tolerance, d);

    Zero = (*Function)(P, beta_R_0, d, Pe, SI);

#if !defined SILENT
    printf("Checking for success in zero finding at %f level of tolerance...\n", Tolerance);
    printf("If a good root has been calculated, then: F(b, d) = 0., within certain tolerance\n");
    printf("F(%7.5f, %7.5f) = %f  +-  %f\n", beta_R_0, d, Zero, Tolerance);
    //Press_Key();
#endif

    return (beta_R_0);
}

float zbrent_Endogen(float (*func)(ParamSet *P, float, float, float, int), ParamSet *P,
		     float x1, float x2, float Pe, int SI, float tol, float NS)
{
	int iter;
	float a=x1,b=x2,c=x2,d,e,min1,min2;
	float fa=(*func)(P, a,NS, Pe,SI),fb=(*func)(P, b,NS, Pe,SI),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		nrerror("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(P, b,NS,Pe,SI);
	}
	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}
#undef ITMAX
#undef EPS





