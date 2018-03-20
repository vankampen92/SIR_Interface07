/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*                      Computing the Analytic Power Spectrum                */
/*                             David Alonso, 2005 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../../SIR.h"
#include "../../../SIR_Analytic_General.h"

int No_of_Points;
double a,d,g,mu,b,Imm;
double d_0,d_1;
double b_0,b_1;
double b_m,b_M;
double Per; 
int POPULATION;
int TIME_DEPENDENT_PARAM;
int NUMBER_of_PLOTS;
int SCAN_WHAT;
int COHERENCE;
double f_semiband;

int main(int argc, char **argv)
{
  double Psi, Fi; /* Populations Fractions */
  double d_Rel, b_Rel;
  int value, point_Damping, point_Stochas, point_Stable;
  double dvalue, coherence, relative_coherence, f_peak;
  int i,j,k;  
  ParamSet Parameters;
  float *Data_1;
  float *Data_2;
  double *x_axis;
  double imm[5] = {0.00001,0.0001,0.001,0.01,0.1};
  char name[22];
  char pre_00[6] = "d_ove_"; 
  char pre_01[6] = "d_coh_"; 
  char pre_10[6] = "b_ove_"; 
  char pre_11[6] = "b_coh_"; 

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  1000; /* Size of the squared matrix to scan: 
			   No_of_Points x No_of_Points */
  NUMBER_of_PLOTS = 5; SCAN_WHAT = 0; TIME_DEPENDENT_PARAM = 0;
  POPULATION = 200000; COHERENCE = 0; f_semiband = 1.E-1; 
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0.; /* Measles as a default */
  b_m = 0.; b_M = 30.; d_0 = 0.; d_1 = 1.;            /* Range of Dimensionless values 
						         to Scan  */
  b_0 = b;  b_1 = 0.25;
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  Per = 365.;
  /* END (Initial settings and default values) * * * * * * * * * * * * * * * * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

  settingParameterStruct(&Parameters);
  modelReport("report.txt");

  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
 
  Data_1 = vector(0,No_of_Points); Data_2 = vector(0,No_of_Points);
  x_axis = dvector(0,No_of_Points);
  
  for(k=0; k<NUMBER_of_PLOTS; k++){
    Parameters.Imm = imm[k]; 
    Parameters.Beta = b; Parameters.Delta = d; Parameters.Gamma = g; 
    Parameters.Alpha = a; Parameters.Mu = mu;  
    /* Making parameters Dimensionless quatities... */ 
    changingTimeScale(&Parameters, 1./g); 
    printf("Dimensionless Parameaters: b/g = %g; g/g = %g; Imm/g = %g\n", 
	   Parameters.Beta, Parameters.Gamma, Parameters.Imm); 
    printf("NOTE: This calculation assumes dimensionless parameters.\n");

    for(i=0; i<No_of_Points; i++){  
      if(SCAN_WHAT == 0){ 
	d_Rel = d_0 + (i+1)* (d_1 - d_0)/(double)No_of_Points;
	if(k==0) x_axis[i] = d_Rel;
	b_Rel = Parameters.Beta;
      }
      else{
	b_Rel = b_m + (i+1)* (b_M - b_m)/(double)No_of_Points;
	if(k==0) x_axis[i] = b_Rel;
	d_Rel = Parameters.Delta;
      }
      re_settingParamStruct(&Parameters, d_Rel, b_Rel);
      
      value = 0;  point_Damping = point_Stochas = point_Stable = 0;
      dvalue = 0.;

      point_Stable = Condition_Stability(&Parameters);
      if(point_Stable == 1){
	point_Damping = 2*Condition(&Parameters, 0.25);
	
	if(point_Damping == 2){
	  //point_Stochas = 4*Exact_Condition_Peak(&Parameters, 1);
	  point_Stochas = 4*Condition(&Parameters, 0.5);
	  dvalue = overall_Amplification(&Parameters, 1);
	}
      }
      /* 
	 The variable $value$ can actually have one out of four values:
	 0: Point is instable;
	 1: Point is stable;
	 3: Point is stable and present damping oscillations;
	 7: Point is stable and present both damping and stochastic amplification;
      */
      value = point_Stable + point_Damping + point_Stochas;
      
      relative_coherence = 0.;
      if(COHERENCE == 1 || COHERENCE == 2){
	if(value == 7){
	  f_peak = Resonance_Frequency_Peak(&Parameters, 1);
	  //coherence = coherence_value(&Parameters, 1, f_peak, f_semiband);
	  //coherence = coherence_value_Simple(f_peak, f_semiband);
	  //relative_coherence = 100. * coherence;   /* /dvalue  */
	  relative_coherence = coherence_value(&Parameters, 1, f_peak, f_semiband);
	}
	Data_2[i] = (float)relative_coherence;
      }
      if(COHERENCE == 0 || COHERENCE == 2){
	Data_1[i] = (float)dvalue;
      }
      if(SCAN_WHAT == 0){	
        printf("d/g = %g\tOverall Amplification: %f\n", x_axis[i], Data_1[i]);
        printf("d/g = %g\tRelative Coherence: %f\n", x_axis[i], Data_2[i]);
      }	
      else{	
        printf("b/g = %g\tOverall Amplification: %f\n", x_axis[i], Data_1[i]);
        printf("b/g = %g\tRelative Coherence: %f\n", x_axis[i], Data_2[i]);
      }
    }
    
    if(SCAN_WHAT == 0){ 
      if(COHERENCE == 0 || COHERENCE == 2){ 
	Saving_to_File_float("dOV_", x_axis, Data_1, No_of_Points, k);
      }
      if(COHERENCE == 1 || COHERENCE == 2){
	Saving_to_File_float("dCO_", x_axis, Data_2, No_of_Points, k);
      }
    }
    else{
      if(COHERENCE == 0 || COHERENCE == 2){ 
	Saving_to_File_float("bOV_", x_axis, Data_1, No_of_Points, k);
      }
      if(COHERENCE == 1 || COHERENCE == 2){
	Saving_to_File_float("bCO_", x_axis, Data_2, No_of_Points, k);
      }
    }  
  }
  
  free_vector(Data_1, 0,No_of_Points); free_dvector(x_axis, 0,No_of_Points);
  free_vector(Data_2, 0,No_of_Points); 
  modelReport("report_End.txt");
  printf("\nEnd of progam\n");
  return (0);
}














