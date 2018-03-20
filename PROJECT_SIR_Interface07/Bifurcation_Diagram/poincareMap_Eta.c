#include "../SIR.h"
#include "../SIR_SEIR.h"

extern double b, b_0, b_1;
extern double b_m, b_M;
extern double Imm;
extern int No_of_Beta;
extern int SEIR, sir;
extern int POPULATION;

extern ParamSIR *Pe;
extern ParamSEIR *Pi;

/* Function involved in the numerical integration of the system */
void poincareMap(int i, int j,
		 float Pop[], int N, float x_i,float x_s, 
		 int No_of_Points, void (*deriva)(double,double[],double[]),
		 double *Time, double **X)
{ 
  FILE *fp;
  int modul, nstep, l, no_Ini, k;
  double Beta, Beta_M; 
  char name[12];
  char file[24];
  double p;
  
  file[0]='\0'; 
  name_Ordered("bifurk.", i, ".", name);
  fitxer(file, name, j, ".dat");      
  fp = fopen(file, "w");

  for(k=0; k<No_of_Beta; k++){
    
    b = b_m + (float)k * (b_M-b_m)/(float)No_of_Beta;
    b_0 = b;  
    Beta_M = 189./365.*b_M + 176./365.*b_M*(1.+b_1); /* Beta Effectiva Maxima */
    Beta = 189./365.*b_0 + 176./365.*b_0*(1.+b_1); /* Beta Effectiva */
    p = (1. - Beta/Beta_M);
    
    if(SEIR == 1){ Pe->b = b; Pe->b_0 = b_0; Pe->Imm = (1. - p) * Imm;}
    if(sir == 1) { Pi->b = b; Pi->b_0 = b_0; Pi->Imm = (1. - p) * Imm;}
    
    integration_Double(Pop, N, x_i, x_s, No_of_Points, deriva, Time, X);
    
    printf(" (1. - p) = %g\n", 1.-p);    
    printf("...........................................%g, %g, %g\n", b_0, b_1, Beta);
    if(SEIR == 1)
      printf("...........................................%g, %g\n", Imm, Pe->Imm);
    if(sir == 1)
      printf("...........................................%g, %g\n", Imm, Pi->Imm);
    
    no_Ini = No_of_Points - 20;
    for (l=no_Ini; l<No_of_Points; l++){
      if(SEIR == 1){
	fprintf(fp, "%g\t%g\n", Beta, X[l][2]/(double)POPULATION);
	if(l > No_of_Points-5) printf("(time: %g, %g)\n",Time[l], X[l][2]/(double)POPULATION);
      }
      else if(sir == 1){
	fprintf(fp, "%g\t%g\n", Beta, X[l][1]/(double)POPULATION);
	if(l > No_of_Points-5) printf("(time: %g, %g)\n",Time[l], X[l][1]/(double)POPULATION);
      }
      else{
	fprintf(fp, "%g\t%g\n", Beta, X[l][1]);
	if(l > No_of_Points-5) printf("(time: %g, %g)\n",Time[l], X[l][1]);
      }
    }
    /*  printf("...........................................%g, %g, %g\n", b_0, b_1, Beta); */
    /*  printf("(%g, %g), (%g, %g), (%g, %g)\n",  */
    /*  Time[No_of_Points-3], X[No_of_Points-3][1],  */
    /*  Time[No_of_Points-2], X[No_of_Points-2][1], */
    /*  Time[No_of_Points-1], X[No_of_Points-1][1]); */
    printf("...\n\n");
  }
  fclose(fp);
}












