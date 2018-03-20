/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*                      Computing the Analytic Power Spectrum                */
/*                                                                           */ 
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
  int value, point_Damping, point_Stochas, point_Stable;
  int i,j,k;  
  ParamSet Parameters;
  FILE *fp, *F_TEX;
  float **Data;
  unsigned char **Cathegory;
  float Max_z, Min_z;
  double imm[5] = {0.00001,0.0001,0.001,0.01,0.1};
  int Label[4] = {2,6,10,14};
  char nameD[22];
  char nameC[22];
  char preC_S[5] = "imS.";
  char preC_I[5] = "imI.";
  char suf[5] = ".eps";

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  1000;
  POPULATION = 100000; SCAN_WHAT = 0; NUMBER_of_PLOTS = 5;
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0.; /* Measles as a default */
  b_0 = b; b_1 = 0.25;
  b_m = 0.; b_M = 30.; d_0 = 0.; d_1 = 1.;            /* Range of Dimensionless values 
						         to Scan                         */
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  /* END (Initial settings and default values) * * * * * * * * * * * * * * * * * * * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);
  
  imm[0] = Imm;
  settingParameterStruct(&Parameters);
  modelReport("report.txt");
  Data = matrix(0,No_of_Points, 0,No_of_Points);
  Cathegory = cmatrix(0,No_of_Points, 0,No_of_Points);
  
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
 
  F_TEX = fopen("./Fig_eps/test_Fig.tex", "w"); 
  printf("Make sure the directory ""./Fig_eps/"" exits!!!\n");
  printf("...otherwise you get a more than beautiful ""Segmentation fault""\n\n");
  preambul_TEX_File(F_TEX);
  
  for(k=0; k<NUMBER_of_PLOTS; k++){
    Parameters.Imm = imm[k]; 
    Parameters.Beta = b; Parameters.Delta = d; Parameters.Gamma = g;  
    Parameters.Alpha = a; Parameters.Mu = mu; 
    /* Making parameters Dimensionless quatities... */ 
    changingTimeScale(&Parameters, 1./g);
    printf("Dimensionless Parameaters: b/g = %g; g/g = %g; Imm/g = %g\n", 
	   Parameters.Beta, Parameters.Gamma, Parameters.Imm); 
    printf("Calculating the borders in the parameter space. \n");
    printf("NOTE: This calculation assumes dimensionless parameters.\n");
    if(k == 0){ 
      if(SCAN_WHAT == 0)
	fp = fopen("region.S.dat", "w"); 
      else
	fp = fopen("region.I.dat", "w");
    }
    for(i=0; i<No_of_Points; i++){
      for(j=0; j<No_of_Points; j++){
	
	value = point_Damping = point_Stochas = point_Stable = 0;
	
	d_Rel = d_0 + (i+1)* (d_1 - d_0)/(double)No_of_Points;
	b_Rel = b_m + (j+1)* (b_M - b_m)/(double)No_of_Points;
	
	re_settingParamStruct(&Parameters, d_Rel, b_Rel);
	//eigen_Values(&Parameters);
	
	point_Stable = Condition_Stability(&Parameters);
	if(point_Stable == 1){
	  point_Damping = 2*Condition(&Parameters, 0.25);
	  
	  if(point_Damping == 2){
	    //point_Stochas = 4*Condition(&Parameters, 0.5);
	    if(SCAN_WHAT == 0){
	      point_Stochas = 4*Exact_Condition_Peak(&Parameters, 0);
	    }
	    else{
	      point_Stochas = 4*Exact_Condition_Peak(&Parameters, 1);
	    }	    
	  }
	}
	value = point_Stable + point_Damping + point_Stochas;
	/*
	  value = 0: Point is instable;
	  value = 1: Point is stable;
	  value = 3: Point is stable and present damping oscillations;
	  value = 7: Point is stable and present both daming and stochastic amplification;
	*/
	switch(value)
	  {
	  case 0: 
	    Cathegory[i][j] = (unsigned char)Label[0];
	    break;  
	  case 1:
	    Cathegory[i][j] = (unsigned char)Label[1];
	    break; 
	  case 3:
	    Cathegory[i][j] = (unsigned char)Label[2];
	    break; 
	  case 7:
	    Cathegory[i][j] = (unsigned char)Label[3];
	    break; 
	  default:
	      printf("  Invalid cathegory: value = %d\n", value);
	      exit(0);
	  }
	Data[i][j] = (float)value;
	if(k == 0) fprintf(fp,"%g\t%g\t%d\n", d_Rel, b_Rel, value);
	//printf("%g\t%g\t%d\n", d_Rel, b_Rel, value);
      }
    }
    
    if(k == 0) fclose(fp);
    
    /* Calculating Scales  and saving eps data */
    vertical_Scale(Data, No_of_Points, &Max_z, &Min_z);
    
    if(SCAN_WHAT == 0){
      name_Ordered(preC_S, k, suf, nameC);
    }
    else{
      name_Ordered(preC_I, k, suf, nameC);
    }
    eps_Cathegory_image(nameC, Cathegory, No_of_Points, No_of_Points);
    append_Fig_TEX(F_TEX, nameC, imm[k]);
  }

  free_matrix(Data, 0,No_of_Points, 0,No_of_Points);
  free_cmatrix(Cathegory, 0,No_of_Points, 0,No_of_Points);
  close_TEX_File(F_TEX);
  modelReport("report_End.txt");
  printf("\nEnd of progam\n");
  return (0);
}














