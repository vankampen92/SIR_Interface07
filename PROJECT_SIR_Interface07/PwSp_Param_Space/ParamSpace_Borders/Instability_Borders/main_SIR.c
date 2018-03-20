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
double d_0,d_1;
double b_0,b_1;
int POPULATION;
double factor;

int main(int argc, char **argv)
{
  double Psi, Fi, Prova_Det, Prova_Dis;  /* Populations Fractions */
  float *d_Rel, *B_0, *B_2, *B_4, *B_2_Higher, *B_4_Higher;
  float x,y, Tr_0, Tr_1;
  double *Traca, *Det;
  int i;  
  ParamSet Parameters;

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  1000;
  POPULATION = 100000;
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0.; /* Measles as a default */
  b_0 = 1.; b_1 = 500.; d_0 = 0.001; d_1 = 1.0;
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  factor = 5.; /* The analytic spectrum is calculated form f=0 to f=factor*f_M, where
		  f_M is the frequency at which the power spectrum peaks */
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

  printf("Closed SIR Assumption: Immigration is negligible, b_e = %g\n", Imm);

  settingParameterStruct(&Parameters);
  modelReport("report.txt");
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
   
  /* Dimensionless parameters... */
  /* Making parameters Dimensionless quatities... */ 
  b *= (1./g); d *= (1./g); Imm *= (1./g); mu *= (1./g); a *= (1./g); g = 1.; 
  settingParameterStruct(&Parameters);
  printf("Calculating the borders in the parameter space. \n");
  printf("NOTE: This calculation assumes dimensionless parameters.\n");
  B_2 = vector(0,No_of_Points);        B_4 = vector(0,No_of_Points); 
  B_2_Higher = vector(0,No_of_Points); B_4_Higher = vector(0,No_of_Points); 
  d_Rel = vector(0,No_of_Points);      B_0 = vector(0,No_of_Points);
    
  for(i=0; i<No_of_Points; i++){
    d_Rel[i] = d_0 + (i+1)* (d_1 - d_0)/(float)No_of_Points;
    printf("Boundaries for Stability (Det > 0)\n");
    B_0[i] = d_Rel[i] + 1.;
    printf("beta_hat (Imm = 0.) = %f\n", B_0[i]);
    Prova_Det = Determinant(&Parameters, B_0[i], d_Rel[i]);
    printf("Det(Beta = %f, Delta = %f, Imm = %g) = %f, only zero if Imm = 0.\n", 
	   B_0[i], d_Rel[i], Parameters.Imm, Prova_Det);
    B_0[i] = border_beta_delta(Determinant, &Parameters, B_0[i], 4., 0.001, d_Rel[i]);
    printf("beta_hat (Imm = %g) = %f\n", Parameters.Imm, B_0[i]); //Press_Key();
    printf("\n");
      
    printf("Boundaries for Stochastic Amplification\n");
    B_2[i] = beta_2(d_Rel[i], -1.0);
    printf("beta_hat (lower root, Imm = 0.) = %f\n", B_2[i]);
    Prova_Dis = Discriminant_2(&Parameters, B_2[i], d_Rel[i]);
    printf("Discriminant(Beta = %f, Delta = %f, Imm = %g) = %f, only zero if Imm = 0.\n",  
	   B_2[i], d_Rel[i], Parameters.Imm, Prova_Dis);
    B_2[i] = border_beta_delta(Discriminant_2, &Parameters, B_0[i], 4., 0.001, d_Rel[i]);
    printf("beta_hat (lower root, Imm = %g) = %f\n", Parameters.Imm, B_2[i]); //Press_Key();

    B_2_Higher[i] = beta_2(d_Rel[i],  1.0);
    printf("beta_hat (higher root, Imm = 0.) = %f\n", B_2_Higher[i]);
    Prova_Dis = Discriminant_2(&Parameters, B_2[i], d_Rel[i]);
    printf("Discriminant(Beta = %f, Delta = %f, Imm = %g) = %f, only zero if Imm = 0.\n", 
	   B_2[i], d_Rel[i], Parameters.Imm, Prova_Dis);
    B_2_Higher[i] = border_beta_delta(Discriminant_2, &Parameters, 4., (float)b_1, 0.001, d_Rel[i]);
    printf("beta_hat (higher root, Imm = %g) = %f\n", Parameters.Imm, B_2_Higher[i]); //Press_Key();

    printf("\n");
    printf("Boundaries for Damped Oscillations\n");
    B_4[i] = beta_4(d_Rel[i], -1.0);
    printf("beta_hat (lower root, Imm = 0.) = %f\n", B_4[i]);
    Prova_Dis = Discriminant_4(&Parameters, B_4[i], d_Rel[i]);
    printf("Discriminant(Beta = %f, Delta = %f, Imm = %g) = %f, only zero if Imm = 0.\n", 
	   B_4[i], d_Rel[i], Parameters.Imm, Prova_Dis);
    B_4[i] = border_beta_delta(Discriminant_4, &Parameters, B_0[i], 4., 0.001, d_Rel[i]);
    printf("beta_hat (lower root, Imm = %g.) = %f\n", Parameters.Imm, B_4[i]); //Press_Key();

    B_4_Higher[i] = beta_4(d_Rel[i],  1.0);
    printf("beta_hat (higher root, Imm = 0.) = %f\n", B_4_Higher[i]);
    Prova_Dis = Discriminant_4(&Parameters, B_4[i], d_Rel[i]);
    printf("Discriminant(Beta = %f, Delta = %f, Imm = %g) = %f, only zero if Imm = 0.\n", 
	   B_4[i], d_Rel[i], Parameters.Imm, Prova_Dis);
    B_4_Higher[i] = border_beta_delta(Discriminant_4, &Parameters, 4., (float)b_1, 0.001, d_Rel[i]);
    printf("beta_hat (higher root, Imm = %g) = %f\n", Parameters.Imm, B_4_Higher[i]); //Press_Key();
    printf("\n\n");
  }
  Saving_to_File_float_float("b_", d_Rel, B_0, No_of_Points, 0);
  Saving_to_File_float_float("b_", d_Rel, B_2, No_of_Points, 2);
  Saving_to_File_float_float("b_", d_Rel, B_4, No_of_Points, 4);
  Saving_to_File_float_float("b_H_", d_Rel, B_2_Higher, No_of_Points, 2);
  Saving_to_File_float_float("b_H_", d_Rel, B_4_Higher, No_of_Points, 4);
  
  free_vector(B_2, 0,No_of_Points);          free_vector(B_4, 0,No_of_Points); 
  free_vector(B_2_Higher, 0,No_of_Points);   free_vector(B_4_Higher, 0,No_of_Points); 
  free_vector(d_Rel, 0,No_of_Points);        free_vector(B_0, 0,No_of_Points);
  
  /* Calculating the curve Det = 1./4 * Tr^2 */
  Traca = dvector(0,No_of_Points); Det = dvector(0,No_of_Points); 
  Tr_1 = 0.1; Tr_0 = -0.1; 
  for(i=0; i<No_of_Points; i++){
      Traca[i] = Tr_0 + (i+1)* (Tr_1 - Tr_0)/(float)No_of_Points;
      Det[i] = 1./2. *Traca[i]*Traca[i];
  }
  Saving_to_File_double("det_", Traca, Det, No_of_Points, 0);
  free_dvector(Traca, 0,No_of_Points);       free_dvector(Det, 0,No_of_Points); 
  
  modelReport("report.txt");
  printf("\nEnd of progam\n");
  return (0);
}














