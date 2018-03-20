/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                       SIR MODEL: Bifurcation Diagram                      */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                            David Alonso, 2000 (c)                         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../SIR.h"
#include "../SIR_SEIR.h"
float x_i, x_s;
float p_1, p_2;
double a;
double mu;
double b;
double b_0, b_1;
double Per;
double d;
double si;
double g;
double Imm;
int MODEL;
int No_of_Points; 
int TIME_DEPENDENT_PARAM;
int POPULATION;
float timeFactor; /* timeFactor is defined as the 1/U, 
		     where U is temporal unit used in days 
		  */
int No_of_Lattice;
double b_m, b_M;
int No_of_Beta;
int SYSTEMATIC;
int SEIR, sir;

ParamSIR *Pi;
ParamSEIR *Pe;

int main(int argc, char **argv)
{
  int N;     /* Number of equations */
  int i,j, Count;
  float *P;
  ParamSet Par;
  ParamSEIR P_SEIR;
  ParamSIR P_SIR;
  FILE *fp_0, *fp_1;
  double f_0, f_1, Beta_eff;
  float p_3,p_4;
  double *Time;
  double **X;

  /* Default values */
  default_values_Measles_02(); 
  if(SEIR == 1)
    values_Measles_02_SEIR(&P_SEIR);
  else if(sir == 1)
    values_Measles_02_SIR(&P_SIR);

  if(argc>1) ArgumentControl(argc,argv);  
  modelReport("report.txt");
  
  if(MODEL == 0) {sir = 0; SEIR = 0; N = 2;}
  else if(MODEL == 1) {sir = 1; SEIR = 0; N = 3; values_Measles_02_SIR(&P_SIR);} 
  else if(MODEL == 2) {sir = 0; SEIR = 1; N = 4; values_Measles_02_SEIR(&P_SEIR);}
  else {printf("Some error in MODEL control variable \n"); exit(0);}
  
  P = vector(0,N-1);
  Time = dvector(0,No_of_Points);   
  X = dmatrix(0,No_of_Points, 0,N-1);
  /* No_of_Lattice * No_of_Lattice bifurcation diagrams are calculated */
  Count = 1;
  for(i=0; i<No_of_Lattice; i++)
    for(j=0; j<No_of_Lattice; j++){
      if(SYSTEMATIC == 0){
	/*  Systematic Exploration of intial conditions     */
	p_1 = (float)(i+1)*0.1/(float)No_of_Lattice;   
	p_2 = (float)(j+1)*0.0001/(float)No_of_Lattice;
	if(SEIR == 1){
	  p_1 = (float)(i+1)*0.1/(float)No_of_Lattice * (float)POPULATION;   
	  p_2 = 0.0001*POPULATION;
	  p_3 = (float)(j+1)*0.0001/(float)No_of_Lattice * (float)POPULATION;
	  p_4 = (float)POPULATION - p_1 - p_2 - p_3; 
	}
	else if(sir == 1){
	  p_1 = (float)(i+1)*0.1/(float)No_of_Lattice * (float)POPULATION;   
	  p_2 = (float)(j+1)*0.0001/(float)No_of_Lattice * (float)POPULATION;
	  p_3 = (float)POPULATION - p_1 - p_2; 
	}
      }
      else{
	/*  Random exploration of initial conditions */
	p_1 = drand48()*0.1;   
	p_2 = drand48()*0.0001;
	if(SEIR == 1){
	  p_1 = drand48() * 0.1 * (float)POPULATION;   
	  p_2 = drand48() * 0.0001 * (float)POPULATION;
	  p_3 = drand48() * 0.0001 * (float)POPULATION;
	  p_4 = (float)POPULATION - p_1 - p_2 - p_3; 
	}
	else if(sir == 1){
	  p_1 = drand48() * 0.1 * (float)POPULATION;   
	  p_2 = drand48() * 0.0001 * (float)POPULATION;
	  p_3 = (float)POPULATION - p_1 - p_2; 
	}
      }
      
      P[0] = p_1;  P[1] = p_2; 
      if(sir == 1) P[2] = p_3;
      if(SEIR == 1){
	P[2] = p_3; P[3] = p_4;
      }
      
      if(MODEL == 0)
	poincareMap(i,j, P, N, x_i, x_s, No_of_Points, deriva, Time, X);
      else if(MODEL == 1)
	poincareMap(i,j, P, N, x_i, x_s, No_of_Points, derivaSIR, Time, X);
      else if(MODEL == 2)
	poincareMap(i,j, P, N, x_i, x_s, No_of_Points, derivaSEIR, Time, X);
      else {printf("Some error in MODEL control variable \n"); exit(0);}
      
      printf(" %dth initial Condition (out of %d): (s=%f, i=%f)\n", 
	     ++Count, No_of_Lattice*No_of_Lattice, p_1,p_2);
    }
  
  printf("Successful Poincare Map!\n");
  free_vector(P, 0, N-1);
  free_dmatrix(X, 0,No_of_Points, 0,N-1);
  free_dvector(Time, 0, No_of_Points);
  return(0);
}
