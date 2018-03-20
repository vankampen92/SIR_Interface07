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
int NUMBER_of_PLOTS;
int SCAN_WHAT;

int main(int argc, char **argv)
{
  double Psi, Fi; /* Populations Fractions */
  double d_Rel, b_Rel;
  int value_S, value_I, point_Damping, point_Stochas_S, point_Stochas_I, point_Stable;
  int i,j,k;  
  ParamSet Parameters;
  float *Data_S, *Data_I, *Data_D, *Data_M;
  double *x_axis;
  double imm[5] = {0.00001,0.0001,0.001,0.01,0.1};
  char name[22];
  float f_peak;

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  1000; /* Size of the squared matrix to scan: 
			   No_of_Points x No_of_Points */
  NUMBER_of_PLOTS = 5; SCAN_WHAT = 0;
  POPULATION = 200000;
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
  Data_S = vector(0,No_of_Points);    Data_I = vector(0, No_of_Points);
  Data_D = vector(0, No_of_Points);   Data_M = vector(0, No_of_Points);    
  x_axis = dvector(0,No_of_Points);       
  /* In the first plot, always there will be imm[0] = Imm, which can 
     be introduced as an input from the command line, for example: 
     >> APWS -e 5.e-4 */
  imm[0] = Imm;
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
 
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
      
      value_S = 0; value_I = 0;
      Data_S[i] = 0.; Data_I[i] = 0.; Data_D[i] = 0.; Data_M[i] = 0.;  
      point_Stable = point_Damping = point_Stochas_S = point_Stochas_I = 0; 
      
      if(SCAN_WHAT == 0){ 
	d_Rel = d_0 + (i+1)* (d_1 - d_0)/(double)No_of_Points;
	x_axis[i] = d_Rel;
	b_Rel = Parameters.Beta;
      }
      else{
	b_Rel = b_m + (i+1)* (b_M - b_m)/(double)No_of_Points;
	x_axis[i] = b_Rel;
	d_Rel = Parameters.Delta;
      }
      
      re_settingParamStruct(&Parameters, d_Rel, b_Rel);
           
      point_Stable = Condition_Stability(&Parameters);
      if(point_Stable == 1){
	point_Damping = 2*Condition(&Parameters, 0.25);
	
	if(point_Damping == 2){
	  point_Stochas_I = 4*Exact_Condition_Peak(&Parameters, 1);
	  point_Stochas_S = 4*Exact_Condition_Peak(&Parameters, 0);
       
	  //point_Stochas = 4*Condition(&Parameters, 0.5);
	  if(Condition(&Parameters, 0.5) == 1){
	    f_peak = (float)Resonance_Frequency(&Parameters)*Per*g;
	    Data_M[i] =  1./f_peak; /* Endegenous period in years */ 
	  }
	  f_peak = (float)Damping_Frequency(&Parameters)*Per*g;
	  Data_D[i] =  1./f_peak; /* Endegenous period in years */ 
	}
      }
      /* 
	 The variable $value$ can actually have one out of four values:
	 0: Point is instable;
	 1: Point is stable;
	 3: Point is stable and present damping oscillations;
	 7: Point is stable and present both damping and stochastic amplification;
      */
      value_S = point_Stable + point_Damping + point_Stochas_S;
      value_I = point_Stable + point_Damping + point_Stochas_I;
      
      if(value_I == 7){
	f_peak = (float)Resonance_Frequency_Peak(&Parameters, 1)*Per*g;
	Data_I[i] =  1./f_peak; /* Endegenous period in years */ 
      }
      if(value_S == 7){
	f_peak = (float)Resonance_Frequency_Peak(&Parameters, 0)*Per*g;
	Data_S[i] =  1./f_peak; /* Endegenous period in years */ 
      }
      if(SCAN_WHAT == 0)
	printf("d/g = %g\tPeriod (years) = %f\n", x_axis[i], Data_I[i]);
      else
	printf("b/g = %g\tPeriod (years) = %f\n", x_axis[i], Data_I[i]);
    } 

    if(SCAN_WHAT == 0){
      Saving_to_File_float("d_EI_", x_axis, Data_I, No_of_Points, k);
      Saving_to_File_float("d_ES_", x_axis, Data_S, No_of_Points, k);
      Saving_to_File_float("d_Da_", x_axis, Data_D, No_of_Points, k);
      Saving_to_File_float("d_Mc_", x_axis, Data_M, No_of_Points, k);
    }
    else{
      Saving_to_File_float("b_EI_", x_axis, Data_I, No_of_Points, k); 
      Saving_to_File_float("b_ES_", x_axis, Data_S, No_of_Points, k); 
      Saving_to_File_float("b_Da_", x_axis, Data_D, No_of_Points, k); 
      Saving_to_File_float("b_Mc_", x_axis, Data_M, No_of_Points, k); 
    } 
  }	
  
  modelReport("report_End.txt");
  printf("\nEnd of progam\n");
  free_vector(Data_S, 0,No_of_Points); free_vector(Data_I, 0, No_of_Points);
  free_dvector(x_axis, 0,No_of_Points);
  free_vector(Data_D, 0, No_of_Points); free_vector(Data_M, 0, No_of_Points);
  return (0);
}














