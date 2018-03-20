/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*                      Computing the Analytic Power Spectrum                */
/*                             David Alonso, 2005 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../SIR.h"
#include "../../SIR_Analytic_General.h"

int No_of_Points;
double a,d,g,mu,b,Imm;
float r;
int POPULATION;
double factor;

int main(int argc, char **argv)
{
  double Psi, Fi, Prova_Det, Prova_Dis;  /* Populations Fractions */
  float *f,*p;
  int i;  
  double Fw, f_r,f_d, fp_S, fp_I, f_semi;
  double P_S, P_I;
  double overall_S, overall_I;
  double overall_S_analytic, overall_I_analytic;
  double coherence_S_analytic, coherence_I_analytic;
  double coherence_S, coherence_I;
  float over,cohe;
  FILE *FP;
  ParamSet Parameters;

  default_values();
  if(argc>1) ArgumentControl(argc,argv);
  settingParameterStruct(&Parameters);
  modelReport("report.txt");
  
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters);
  printf("\n");
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);

  if(Condition(&Parameters, 0.25) == 1){
    f_r = Resonance_Frequency(&Parameters);
    f_d = Damping_Frequency(&Parameters);
    fp_S = Resonance_Frequency_Peak(&Parameters, 0);
    fp_I = Resonance_Frequency_Peak(&Parameters, 1);  
    
    P_S = Maximum_Power(&Parameters,0);
    P_I = Maximum_Power(&Parameters,1);
    
    printf("\nPeak Frequency (S): fp_S = %g\tPeak Frequency (I): fp_I = %g\n",
	   fp_S, fp_I);
    fp_S = Resonance_Frequency_Peak_Approx(&Parameters, 0);
    fp_I = Resonance_Frequency_Peak_Approx(&Parameters, 1);  
    printf("\nPeak Frequency Aprox.(S): fp_S = %g\tPeak Frequency Approx.(I): fp_I = %g\n",
	   fp_S, fp_I);
    printf("\nResonance Frequency: f_r = %g\tDamping Frequency: f_d = %g\n",
	   f_r, f_d);
    printf("\nMaximum Spectrum Power (S): %g\tMaximum Spectrum Power (I): %g\n",
	   P_S, P_I);
    
  /* Saving Power Spectral Peaks... */
  /* True Peak */
  fp_S = Resonance_Frequency_Peak(&Parameters, 0);
  fp_I = Resonance_Frequency_Peak(&Parameters, 1);
  FP = fopen("PeakFr_S.dat", "w");
  fprintf(FP, "%g\t%g\n", fp_S, 0.); fprintf(FP, "%g\t%g\n", fp_S, P_S);
  fclose(FP);
  FP = fopen("PeakFr_I.dat", "w");
  fprintf(FP, "%g\t%g\n", fp_I,0.); fprintf(FP, "%g\t%g\n", fp_I, P_I);
  fclose(FP);
  /* Damping Peaks */
  FP = fopen("PeakDamp_S.dat", "w");
  fprintf(FP, "%g\t%g\n", f_d, 0.); fprintf(FP, "%g\t%g\n", f_d, P_S);
  fclose(FP);
  FP = fopen("PeakDamp_I.dat", "w");
  fprintf(FP, "%g\t%g\n", f_d,0.); fprintf(FP, "%g\t%g\n", f_d, P_I);
  fclose(FP);
  /* Resonance Peaks */
  FP = fopen("PeakReson_S.dat", "w");
  fprintf(FP, "%g\t%g\n", f_r, 0.); fprintf(FP, "%g\t%g\n", f_r, P_S);
  fclose(FP);
  FP = fopen("PeakReson_I.dat", "w");
  fprintf(FP, "%g\t%g\n", f_r,0.); fprintf(FP, "%g\t%g\n", f_r, P_I);
  fclose(FP);
  /* End saving Peaks */
  }
  else{
    printf("\n  The condition Tra^2(A) - 4.*Det(A) < 0. is not fulfilled.\n");
    printf("  The power spectra do not present a meaningful peak.\n");
    if(Condition_Stability(&Parameters) == 1){    
      printf("  The fix point is stable but is not an attractive spiral.\n");
    }
    else{
      printf("  The fix point is not even stable\n");
    }
    f_r = 0.005;
    factor = 1.;
  }
  
  /* Saving Power Spectral Density for Susceptibles */
  f = vector(0,No_of_Points); p = vector(0,No_of_Points); 
  for(i=0; i<No_of_Points; i++){
    f[i] = factor * f_r/(double)No_of_Points * (double)i;
    Fw   = 2.*M_PI*f[i];
    p[i] = power_Spectrum(&Parameters, 0, Fw);
  }
  Saving_to_File_float_float("pwsp", f, p, No_of_Points, 0, 0);
  
  /* Saving Power Spectral Density for Infectives */
  f = vector(0,No_of_Points); p = vector(0,No_of_Points); 
  for(i=0; i<No_of_Points; i++){
    f[i] = factor * f_r/(double)No_of_Points * (double)i;
    Fw   = 2.*M_PI*f[i];
    p[i] = power_Spectrum(&Parameters, 1, Fw);
  }
  Saving_to_File_float_float("pwsp", f, p, No_of_Points, 1, 0);
  
  /* Calculating the area under the spectral density curve
     Integratring spectral densities overall all frequencies */
  printf("\n Area under the spectral density curve:\n");

  overall_S_analytic = overall_Analytic(&Parameters, 0);
  overall_I_analytic = overall_Analytic(&Parameters, 1);
  printf(" Susceptible (Analytic): %g\n", overall_S_analytic);
  printf(" Infective   (Analytic): %g\n", overall_I_analytic);

  overall_S =  overall_Amplification(&Parameters, 0);
  overall_I =  overall_Amplification(&Parameters, 1);
  printf(" Susceptible (Numeric): %g\n", overall_S);
  printf(" Infective   (Numeric): %g\n", overall_I);

  /* Calculating the area around the spectral density peak
     Integratring spectral densities around the peak */
  printf("\n Coherence of the spectral density curve:\n");
  
  f_semi = 2. * (double)r * fp_S; over = overall_Analytic(&Parameters, 0);
  coherence_S_analytic = 100.*coherence_Analytic(&Parameters, 0, f_semi, fp_S)/over;
  f_semi = 2. * (double)r * fp_I; over = overall_Analytic(&Parameters, 1);
  coherence_I_analytic = 100.*coherence_Analytic(&Parameters, 1, f_semi, fp_I)/over;
  printf(" Susceptible (Analytic): %g\n", coherence_S_analytic);
  printf(" Infective   (Analytic): %g\n", coherence_I_analytic);

  f_semi = 2. * (double)r * fp_S;
  coherence_S =  coherence_semi_analytic(&Parameters, 0, fp_S, f_semi);
  f_semi = 2. * (double)r * fp_I;
  coherence_I =  coherence_semi_analytic(&Parameters, 1, fp_I, f_semi);
  printf(" Susceptible (Numeric): %g\n", coherence_S);
  printf(" Infective   (Numeric): %g\n", coherence_I);
  
  
  printf("\n Overall amplification and coherence of the spectral density curve (2nd. cal):\n");
  f_semi = 2. * (double)r * fp_S;
  printf(" Power Spectral density for periods (in years) between %g and %g\n", 
	 1./365./(fp_S-0.5*f_semi), 1./365./(fp_S+0.5*f_semi));
  overall_Amplification_Coherence(&Parameters, 0, fp_S, f_semi,
				  &over, &cohe);
  printf(" Susceptible (Numeric): Amplifiation: %g\tCoherence: %g\n", over,cohe);

  f_semi = 2. * (double)r * fp_I;
  printf(" Power Spectral density for periods (in years) between %g and %g\n", 
	 1./365./(fp_I+0.5*f_semi), 1./365./(fp_I-0.5*f_semi));
  overall_Amplification_Coherence(&Parameters, 1, fp_I, f_semi,
				  &over, &cohe);
  printf(" Infective (Numeric): Amplifiation: %g\tCoherence: %g\n", over,cohe);

  free_vector(f, 0,No_of_Points); free_vector(p, 0,No_of_Points);   
  printf("\n");
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}













