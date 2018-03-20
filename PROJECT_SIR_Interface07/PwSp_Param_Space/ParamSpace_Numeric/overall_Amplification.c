/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*                                                                           */
/*                      Computing the Analytic Power Spectrum                */ 
/*                         and the overall amplification                     */
/*                            (mean squared deviation)                       */
/*                                                                           */
/*                             David Alonso, 2005 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../SIR.h"
#include "../../SIR_Analytic_General.h"

float power_Spec_Func(float Fw);
float power_Spec_Func_overall(float Fw);
ParamSet *Parameters;
int SI_Func;
double overall;

double overall_Amplification(ParamSet *P, int SI)
{
  double Psi, Fi, Prova_Det, Prova_Dis;  /* Populations Fractions */
  double factor, f_r;
  
  float AINF, spcden_AINF, EPS;
  static int Count = 0;

  Parameters = P;
  SI_Func = SI;
  
  if(Condition(Parameters, 0.25) != 1){
    printf("  Stochastic amplification can only potentially exist and be calculated\n");
    printf("  wherever there are damped oscillations to a stable point\n");
    printf("\n  The condition Tra^2(A) - 4.*Det(A) < 0. is not fulfilled.\n");
    printf("  The power spectra do not present a meaningful peak.\n");
    if(Condition_Stability(Parameters) == 1){    
      printf("  The fix point is stable but is not an attractive spiral.\n");
    }
    else{
      printf("  The fix point is not even stable\n");
    }
    f_r = 0.005;
    factor = 1.;
    overall = 0.;
  }
  else{
    /* Maximum freq is A cicles per time unit: AINF = A * 2. * M_PI; */ 
    AINF =  13.*2.*M_PI;  spcden_AINF = power_Spec_Func(AINF);
    EPS = 1.0E-6;
   
    overall = qromb_Accuracy(power_Spec_Func, 0., AINF, EPS);
    //overall = qtrap_EPS(power_Spec_Func, 0., AINF, EPS);
    
    //printf("  ... ... Call #%d: Total__Power(b/g = %g, d/g %g) = %g\tSpcDen(w=AINF)=%f\n", 
    //	   ++Count, Parameters->Beta, Parameters->Delta, overall, spcden_AINF);
  }
  return (overall);
}

double coherence_value(ParamSet *P, int SI, double f_peak, double f_semi)
{
  float ASUP, AINF, EPS, spcden_peak;
  static int Count = 0;
  double coherence;
  
  Parameters = P;
  SI_Func = SI;

  ASUP = (f_peak + 0.5*f_semi)*2.*M_PI;
  AINF = (f_peak - 0.5*f_semi)*2.*M_PI;
  EPS  = 5.0E-7;

  if(overall > 0.) {
    //coherence = qtrap_EPS(power_Spec_Func, AINF, ASUP, EPS);
    coherence = qromb_Accuracy(power_Spec_Func_overall, AINF, ASUP, EPS);
    if(coherence > 100.){ 
      printf("Some error in the estimation of the the amplification and coherence\n");
      coherence = 0.;
    }
    //spcden_peak = power_Spec_Func(f_peak);
    //printf("  ... ... Call #%d: Coherence__Power(b/g = %g, d/g %g) = %g\tSpcDen(w=peak)=%f\n", 
    //	 ++Count, Parameters->Beta, Parameters->Delta, coherence, spcden_peak);
  }
  else{
    coherence = 0.;
  }

  return(coherence);
}

double coherence_semi_analytic(ParamSet *P, int SI, double f_peak, double f_semi)
{
  float ASUP, AINF, EPS, spcden_peak;
  static int Count = 0;
  double coherence;
  
  Parameters = P;
  SI_Func = SI;

  ASUP = (f_peak + 0.5*f_semi)*2.*M_PI;
  AINF = (f_peak - 0.5*f_semi)*2.*M_PI;
  EPS  = 5.0E-7;
  
  overall = overall_Analytic(P,SI);
  if(overall > 0.) {   
    //coherence = qtrap_EPS(power_Spec_Func, AINF, ASUP, EPS);
    coherence = qromb_Accuracy(power_Spec_Func_overall, AINF, ASUP, EPS);
    if(coherence > 100.){ 
      printf("Some error in the estimation of the the amplification and coherence\n");
      coherence = 0.;
    }
    spcden_peak = power_Spec_Func(f_peak);
    printf("  ... ... Call #%d: Coherence__Power(b/g = %g, d/g %g) = %g\tSpcDen(w=peak)=%f\n", 
	   ++Count, Parameters->Beta, Parameters->Delta, coherence, spcden_peak);
  }
  else{
    coherence = 0.;
  }

  return(coherence);
}

void overall_Amplification_Coherence(ParamSet *P, int SI, double f_peak, double f_semi,
				     float *over, float *cohe)
{
  double Psi, Fi, Prova_Det, Prova_Dis;  /* Populations Fractions */
  double factor, f_r; 
  float AINF, ASUP, EPS;
 
  Parameters = P;
  SI_Func = SI;
  
  if(Condition(Parameters, 0.25) != 1){
    printf("  Stochastic amplification can only potentially exist and be calculated\n");
    printf("  wherever there are damped oscillations to a stable point\n");
    printf("\n  The condition Tra^2(A) - 4.*Det(A) < 0. is not fulfilled.\n");
    printf("  The power spectra do not present a meaningful peak.\n");
    if(Condition_Stability(Parameters) == 1){    
      printf("  The fix point is stable but is not an attractive spiral.\n");
    }
    else{
      printf("  The fix point is not even stable\n");
    }
    f_r = 0.005;
    factor = 1.;
  }
  else{
    /* Maximum freq is A cicles per time unit: AINF = A * 2. * M_PI; */ 
    AINF =  13.*2.*M_PI;  
    EPS = 1.0E-6;
    (*over) = qromb_Accuracy(power_Spec_Func, 0., AINF, EPS);
    overall = (*over);
    
    if(overall > 0.){
      ASUP = (f_peak + 0.5*f_semi)*2.*M_PI;
      AINF = (f_peak - 0.5*f_semi)*2.*M_PI;
      EPS  = 5.0E-7;
      (*cohe) = qromb_Accuracy(power_Spec_Func_overall, AINF, ASUP, EPS);
      //assert( *cohe <= 100. );
      if( *cohe > 100.) *cohe = 0.;
    }
    else{
      (*cohe) = 0.;
    }
  }
}

float power_Spec_Func(float Fw)
{
  double Fww, Val;
  float value;  
  
  Fww = (double)Fw;
  Val = power_Spectrum(Parameters, SI_Func, Fww);
  
  value = (float)Val;
  
  return value;
}

float power_Spec_Func_overall(float Fw)
{
  double Fww, Val;
  float value;  
  
  Fww = (double)Fw;
  Val = 100. * power_Spectrum(Parameters, SI_Func, Fww) / overall;
  
  value = (float)Val;
  
  return value;
}





























