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
int POPULATION;
float TOLERANCE;
int NUMBER_of_PLOTS;
int No_of_PERIOD_VALUES;
int SCAN_WHAT;

float constrained_Beta(ParamSet *P, 
		       float Beta_0, float Delta, float Eps, int SI, int signe);

int main(int argc, char **argv)
{
  double Psi, Fi, Prova_Det, Prova_Dis;  /* Populations Fractions */
  float *d_Rel, *B_Rel, *B_Inf, *B_Sup;
  float x,y, Tr_0, Tr_1;
  int i,j,k;  
  ParamSet Parameters;
  float Period[100],Per_0;
  float imm[5] = {0.00001,0.0001,0.001,0.01,0.1};
  char name[22];
  float F0, F1, b_inf, b_sup;
  int Minor_Ticks, Major_Ticks;
  float Major_Spacing;
  
  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  1000; 
  No_of_PERIOD_VALUES = 12; NUMBER_of_PLOTS = 1; TOLERANCE = 0.001;
  POPULATION = 100000;
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0.; /* Measles as a default */
  b_0 = 0.; b_1 = 30.; d_0 = 0.; d_1 = 1.0;
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  Per_0 = 0.5; Minor_Ticks = 3; Major_Ticks = 7; Major_Spacing = 1.; 
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

  settingParameterStruct(&Parameters);
  modelReport("report.txt");
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
  
  /* Initiating vector of Periods */
  k=0;
  for(i=0; i<Major_Ticks; i++){
    for(j=0; j<=Minor_Ticks; j++){ 
      Period[k++] = Per_0 + (float)i*Major_Spacing + (float)j*Major_Spacing/((float)Minor_Ticks + 1.);
      printf("Period... %f\n", Period[k-1]);
    }
  }
  Press_Key();
  
  No_of_PERIOD_VALUES = k;
  B_Rel = vector(0,No_of_Points); d_Rel = vector(0,No_of_Points);    
  B_Inf = vector(0,No_of_Points); B_Sup = vector(0,No_of_Points);
  
  for(k=0; k<NUMBER_of_PLOTS; k++){
    Parameters.Imm = imm[k]; 
    Parameters.Beta = b; Parameters.Delta = d; Parameters.Gamma = g; 
    Parameters.Alpha = a; Parameters.Mu = mu;  
    /* Making parameters Dimensionless quatities... */ 
    changingTimeScale(&Parameters, 1./g); 
    printf("Dimensionless Parameaters: b/g = %g; g/g = %g; Imm/g = %g\n", 
	   Parameters.Beta, Parameters.Gamma, Parameters.Imm); 
    printf("NOTE: This calculation assumes dimensionless parameters.\n");
    for(j=0; j<No_of_PERIOD_VALUES; j++){ 
      for(i=0; i<No_of_Points; i++){
	printf("Boundaries for Endogenous Period\n");
	d_Rel[i] = d_0 + (i+1)*(d_1 - d_0)/(float)No_of_Points; /* This is a guess */
	b_inf = constrained_Beta(&Parameters, 
				 (float)b_0, d_Rel[i], TOLERANCE, SCAN_WHAT, 1); 
	b_sup = constrained_Beta(&Parameters, 
				 (float)b_1, d_Rel[i], TOLERANCE, SCAN_WHAT,-1); 
	
	F0 = endogen_Period(&Parameters, b_inf,d_Rel[i], Period[j],SCAN_WHAT); 
	F1 = endogen_Period(&Parameters, b_sup,d_Rel[i], Period[j],SCAN_WHAT);	
	if(F0*F1 > 0){ 
	  printf("Searching for %f > b > %f,\n", b_inf, b_sup);
	  printf("Other values: a = %g, d = %g, mu = %g, g = %g, Imm = %g\n",
		 Parameters.Alpha, Parameters.Delta, Parameters.Mu, Parameters.Gamma, 
		 Parameters.Imm);
	  printf("No value found... Failure bracking the root\n");
	  printf("Function(b_inf=%f,d=%f) = %f\tFunction(b_sup=%f,d=%f) = %f\n", 
		 b_inf,d_Rel[i],F0, b_sup,d_Rel[i],F1);  
	  B_Rel[i] = 0.;
	}
	else
	  B_Rel[i] = border_Period(endogen_Period, &Parameters, b_inf, b_sup,
				   Period[j], SCAN_WHAT, TOLERANCE, d_Rel[i]);
	
	printf("Searching for %f > b > %f,\n", b_inf, b_sup);
	B_Inf[i] = b_inf; B_Sup[i] = b_sup;
	printf("Period(years)= %f\t Delta= %f\t Beta= %f\n", Period[j], d_Rel[i], B_Rel[i]);
      }
      name_Ordered("Contour_", k, "_", name);
      Saving_to_File_float_float(name, d_Rel, B_Rel, No_of_Points, j, 1);
    }
    name_Ordered("B_Inf_", k, "_",name);
    Saving_to_File_float_float(name, d_Rel, B_Inf, No_of_Points, 0, 1);
    name_Ordered("B_Sup_", k, "_",name);
    Saving_to_File_float_float(name, d_Rel, B_Sup, No_of_Points, 0, 1);
  }
  
  free_vector(B_Rel, 0,No_of_Points);    free_vector(d_Rel, 0,No_of_Points);       
  free_vector(B_Inf, 0,No_of_Points);    free_vector(B_Sup, 0,No_of_Points);  
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}

float constrained_Beta(ParamSet *P, 
		       float Beta_0, float Delta, float Eps, int SI, int signe)
{
  /* Beta values should be constrained in the region where there is existence of
     an "interior peak" in the spectral density */  
  int i;
  float Beta;
  
  Beta = Beta_0; 
  P->Beta = Beta_0;
  P->Delta = Delta;
  
  while(Exact_Condition_Peak(P,SI) == 0) P->Beta = P->Beta + signe*Eps;
   
  Beta = P->Beta;
  return(Beta);
}












