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

int main(int argc, char **argv)
{
  double Psi, Fi; /* Populations Fractions */
  double d_Rel, b_Rel;
  int value, point_Damping, point_Stochas, point_Stable;
  int i,j,k;  
  ParamSet Parameters;
  FILE *fp, *F_TEX;
  float **Data;
  float f_peak_S, f_peak_D; 
  float Max_z, Min_z;
  double imm[5] = {0.00001,0.0001,0.001,0.01,0.1};
  char name[22];
  char pre[5] = "imm_";
  char su[5] = ".eps";

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  1000;
  POPULATION = 100000;
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0.; /* Measles as a default */
  b_m = 0.; b_M = 30.; d_0 = 0.; d_1 = 1.;            /* Range of Dimensionless values 
						         to Scan */
  b_0 = b;  b_1 = 0.25;
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  Per = 365.;
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);
  imm[0] = Imm; 
  settingParameterStruct(&Parameters);
  modelReport("report.txt");
  Data = matrix(0,No_of_Points, 0,No_of_Points);
  
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
 
  F_TEX = fopen("./Fig_eps/test_Fig.tex", "w");
  preambul_TEX_File(F_TEX);

  for(k=0; k<1; k++){
    Parameters.Imm = imm[k]; 
    Parameters.Beta = b; Parameters.Delta = d; Parameters.Gamma = g; 
    Parameters.Mu = mu; Parameters.Alpha = a; 
    /* Making parameters Dimensionless quatities... */ 
    changingTimeScale(&Parameters, 1./g);
    printf("Dimensionless Parameaters: b/g = %g; g/g = %g; Imm/g = %g\n", 
	   Parameters.Beta, Parameters.Gamma, Parameters.Imm); 
    printf("Calculating the borders in the parameter space. \n");
    printf("NOTE: This calculation assumes dimensionless parameters.\n");
    
    fp = fopen("paramScan.dat", "w"); 
    for(i=0; i<No_of_Points; i++){
      for(j=0; j<No_of_Points; j++){
	
	point_Stable = point_Stochas = point_Damping = 0;
	
	d_Rel = d_0 + (i+1)* (d_1 - d_0)/(double)No_of_Points;
	b_Rel = b_m + (j+1)* (b_M - b_m)/(double)No_of_Points;
	
	re_settingParamStruct(&Parameters, d_Rel, b_Rel);
	//eigen_Values(&Parameters);
	
	value = 0;
	point_Stable = Condition_Stability(&Parameters);
	if(point_Stable == 1){
	  point_Damping = 2*Condition(&Parameters, 0.25);	
	  if(point_Damping == 2){
	    point_Stochas = 4*Exact_Condition_Peak(&Parameters, 1);
	    //point_Stochas = 4*Condition(&Parameters, 0.5);
	  }
	}
	value = point_Stochas + point_Stable + point_Damping;
	
	if(value == 7){
	  /* There is a distinctive peak from zero in the spectra of 
	     infectives. Data stores the endogenous period in year units. 
	  */
	  f_peak_S = (float)Resonance_Frequency_Peak(&Parameters, 1)*Per*g;
	  /* Period associated to the deterministic damped oscillations
	     to the stable point in year units 
	  */
	  f_peak_D = (float)Imaginary_Eigen_Value(&Parameters)*Per*g/2./M_PI; 	
	  //Data[i][j] = 1./f_peak_D - 1./f_peak_S;
	  Data[i][j] = 1./f_peak_S;
 	}
	else
	  Data[i][j] = 0.;

	fprintf(fp,"%g\t%g\t%f\n", d_Rel, b_Rel, Data[i][j]);
	printf("%g\t%g\t%f\n", d_Rel, b_Rel, Data[i][j]);
      }
    }
    fclose(fp);

    name_Ordered(pre, k, su, name);
    /* Calculating Scales  and saving eps data */
    vertical_Scale(Data, No_of_Points, &Max_z, &Min_z);
    eps_WGW_image(name, Data, No_of_Points, No_of_Points, Max_z, Min_z);

    append_Fig_TEX(F_TEX, name, imm[k]);
  }

  close_TEX_File(F_TEX);
  modelReport("report_End.txt");
  printf("\nEnd of progam\n");
  free_matrix(Data, 0,No_of_Points, 0,No_of_Points);
  
  return (0);
}














