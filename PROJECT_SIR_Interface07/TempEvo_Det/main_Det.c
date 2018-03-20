/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                       SIR MODEL: Numerical Integration                    */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                            David Alonso, 2000 (c)                         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../SIR.h"

float x_i, x_s;
float p_1, p_2;
double a;
double mu;
double b;
double b_0, b_1;
double Per;
double g;
double d;
double Imm;
int N;     /* Number of equations */
int No_of_Points; 
int TIME_DEPENDENT_PARAM;
int POPULATION;
float timeFactor;

int main(int argc, char **argv)
{
  int i,j;
  float *P;
  ParamSet Par;
  FILE *fp_0, *fp_1;
  double f_0, f_1;
  double *Time;
  double **X;

  /* Default values */
  N = 2; /* Number of Equations */
  POPULATION = 100000;
  No_of_Points = 1000; TIME_DEPENDENT_PARAM = 0;
  x_i=0.; x_s=3000.;
  p_1=0.1; p_2=0.0001;
  b = 1.175; a = 0.; d = 5.5e-05; g = 0.077; mu = 0.; Imm = 1.e-05;
  b_0 = b; b_1 = 0.04; timeFactor=1.; Per = 365./timeFactor;
 
  if(argc>1) ArgumentControl(argc,argv);
  
  modelReport("report.txt");

  P = vector(0,N-1);
  Time = dvector(0,No_of_Points);   
  X = dmatrix(0,No_of_Points, 0,N-1);

  P[0] = p_1; P[1] = p_2; b_0 = b;
  Par.Beta = b; Par.Alpha= a; Par.Delta= d; Par.Gamma= g; Par.Mu = mu;
  Par.Imm = Imm;
  /* integration(...) is a function saves "Number_of_Points" of the temporal 
     evolution in a file mfTempEvo.dat */
  
  integration_01(P, N, x_i, x_s, No_of_Points, deriva, Time, X);
  //integration(P, N, x_i, x_s, No_of_Points, deriva);  
  printf("Successful numerical integration!\n\n");

  printf("Writing temporal evolution in population numbers\n");
  fp_1 = fopen("mfTempEvo_No.dat", "w");
  for(i=0; i<No_of_Points; i++){
    p_1 = X[i][0] * POPULATION; p_2 = X[i][1]* POPULATION;
    fprintf(fp_1, "%f\t%f\t%f\n", Time[i], p_1, p_2);
  }
  fclose(fp_1);
  
  
  printf("Stability analysis of the immigration-closed seasonal non-forced system\n");
  Fixed_Points(&Par, &f_0, &f_1); Stability(&Par);

  free_vector(P, 0, N-1);
  free_dmatrix(X, 0,No_of_Points, 0,N-1);
  free_dvector(Time, 0, No_of_Points);
  
  return(0);
}















