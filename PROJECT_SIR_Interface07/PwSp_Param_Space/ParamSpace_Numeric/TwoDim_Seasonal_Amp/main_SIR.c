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
int STO_DET;
int MSD; 
int NUMBER_of_PLOTS;
float timeFactor; /* timeFactor is defined as the 1/U, 
		     where U is temporal unit used in days 
		  */
int main(int argc, char **argv)
{
  double Psi, Fi; /* Populations Fractions */
  double d_Rel, b_Rel;
  int value, point_Damping, point_Stochas, point_Stable;
  int i,j,k;  
  ParamSet Parameters;
  FILE *fp, *F_TEX;
  float **Data;
  float Max_z, Min_z;
  double imm[5] = {0.00001,0.0001,0.001,0.01,0.1};
  char name[22];
  char pre[5] = "imm_";
  char suf[5] = ".eps";

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  1000; /* Size of the squared matrix to scan: 
			   No_of_Points x No_of_Points */
  TIME_DEPENDENT_PARAM = 3;
  POPULATION = 500000; STO_DET = 0; MSD = 1; NUMBER_of_PLOTS = 1;
  timeFactor = 1/13.;
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
  Data = matrix(0,No_of_Points, 0,No_of_Points);
  
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
  
  F_TEX = fopen("./Fig_eps/test_Fig.tex", "w");
  preambul_TEX_File(F_TEX);
 
  for(k=0; k<NUMBER_of_PLOTS; k++){
    Parameters.Imm = imm[k]; g = 1./13.; Per = 365.*timeFactor;
    Parameters.Beta = b; Parameters.Delta = d; Parameters.Gamma = g; 
    Parameters.Alpha = a; Parameters.Mu = mu;  
    /* Making parameters Dimensionless quatities... */ 
    changingTimeScale(&Parameters, 1./timeFactor); 
    printf("Dimensionless Parameaters: b/g = %g; g/g = %g; Imm/g = %g\n", 
	   Parameters.Beta, Parameters.Gamma, Parameters.Imm); 
    printf("NOTE: This calculation assumes dimensionless parameters.\n");

    if(k==0) fp = fopen("paramScan.dat", "w");

    for(i=0; i<No_of_Points; i++){
      for(j=0; j<No_of_Points; j++){
	
	point_Damping = point_Stochas = point_Stable = 0;
	
	d_Rel = d_0 + (i+1)* (d_1 - d_0)/(double)No_of_Points;
	b_Rel = b_m + (j+1)* (b_M - b_m)/(double)No_of_Points;
	
	re_settingParamStruct(&Parameters, d_Rel, b_Rel);
	//eigen_Values(&Parameters);
	
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

	if(value == 3 || value == 7)
	  if(STO_DET == 0)  {Data[i][j] = sto_season_amplif(MSD, &Parameters, t_f, 5, 1);} 
	  else              {Data[i][j] = det_season_amplif(MSD, &Parameters, t_f, 5, 1);}
	  /* There are at least damping oscillations */
 	else
	  Data[i][j] = 0.;
	if(k==0) fprintf(fp,"%g\t%g\t%f\n", d_Rel, b_Rel, Data[i][j]);
      }
    }
    if(k==0) fclose(fp);
    
    name_Ordered(pre, k, suf, name); /* Files are labelled in order: pre_k.suf */
    /* Calculating Scales  and saving eps data */
    vertical_Scale(Data, No_of_Points, &Max_z, &Min_z);
    eps_WGW_image(name, Data, No_of_Points, No_of_Points, Max_z, Min_z);

    append_Fig_TEX(F_TEX, name, imm[k]);
    //ps_image("overall.eps", Data, No_of_Points, Max_z, Min_z);
  }

  close_TEX_File(F_TEX);
  modelReport("report_End.txt");
  printf("\nEnd of progam\n");
  free_matrix(Data, 0,No_of_Points, 0,No_of_Points);
  return (0);
}














