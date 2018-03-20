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
float t_f;      /* Final integration time */
double Per; 
int POPULATION;
int TIME_DEPENDENT_PARAM;
int NUMBER_of_PLOTS;
int SCAN_WHAT;
int STO_DET;
int MSD; 

float timeFactor; /* timeFactor is defined as the 1/U, 
		     where U is temporal unit used in days float
		  */

int main(int argc, char **argv)
{
  double Psi, Fi; /* Populations Fractions */
  double d_Rel, b_Rel;
  int value, point_Damping, point_Stochas, point_Stable;
  int i,j,k;  
  ParamSet Parameters;
  float *Data;
  double *x_axis;
  double imm[5] = {0.00001,0.0001,0.001,0.01,0.1};
  char name[22];
  char pre_0[6] = "d_imm_";
  char pre_1[6] = "b_imm_";

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * * */
  No_of_Points =  1000; /* Size of the squared matrix to scan: 
			   No_of_Points x No_of_Points */
  NUMBER_of_PLOTS = 5; SCAN_WHAT = 0; STO_DET = 0; TIME_DEPENDENT_PARAM = 0;
  POPULATION = 200000; MSD = 0;
  time_Factor = 1./13.;
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0.; /* Measles as a default */
  b_m = 0.; b_M = 30.; d_0 = 0.; d_1 = 1.;            /* Range of Dimensionless values 
						         to Scan  */
  b_0 = b;  b_1 = 0.25;
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  Per = 365.;
  t_f = 500.;  /* t_f in day units (or in dimensionalless units, BE CAREFUL! */
  /* END (Initial settings and default values) * * * * * * * * * * * * * * * * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

  settingParameterStruct(&Parameters);
  modelReport("report.txt");
  Data = vector(0,No_of_Points);
  x_axis = dvector(0,No_of_Points);

  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
 
  for(k=0; k<NUMBER_of_PLOTS; k++){
    Parameters.Imm = imm[k]; g = 1./13.; Per = 365.*time_Factor; /* Period in dimensionless units */
    Parameters.Beta = b; Parameters.Delta = d; Parameters.Gamma = g; 
    Parameters.Alpha = a; Parameters.Mu = mu;  
    /* Making parameters Dimensionless quatities... */ 
    changingTimeScale(&Parameters, 1./time_Factor); 
    
    printf("NOTE: This calculation assumes dimensionless parameters.\n");
    printf("Immigration level = %g\n", Parameters.Imm); 

    for(i=0; i<No_of_Points; i++){
      
      point_Damping = point_Stochas = point_Stable = 0;
      
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
	  point_Stochas = 4*Exact_Condition_Peak(&Parameters, 1);
	  //point_Stochas = 4*Condition(&Parameters, 0.5);
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
      
      if(value == 3 || value == 7){
	/* There are at least deterministic damping oscillations 
	   The five largest peaks are averaged */
	if(STO_DET == 0)  {Data[i] = sto_season_amplif(MSD, &Parameters, t_f, 5, 1);} 
        else              {Data[i] = det_season_amplif(MSD, &Parameters, t_f, 5, 1);}
      }
      else
	Data[i] = 0.;
    }
    
    if(SCAN_WHAT == 0) 
      Saving_to_File_float(pre_0, x_axis, Data, No_of_Points, k);
    else
      Saving_to_File_float(pre_1, x_axis, Data, No_of_Points, k); 
  }
  
  modelReport("report_End.txt");
  printf("\nEnd of progam\n");
  free_vector(Data, 0,No_of_Points); free_dvector(x_axis, 0,No_of_Points);
  return (0);
}














