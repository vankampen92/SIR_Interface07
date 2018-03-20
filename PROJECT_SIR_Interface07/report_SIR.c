#include "SIR.h" 
extern int TIMES;     
extern float STEP_SIZE; 
extern float EPSILON;
extern int Realizations;
extern int Simlength;
extern double o1;
extern double g1;
extern double a1;
extern double mu;
extern double bh1;
extern double bw1;
extern double r1;
extern double f;  /* Water Infectivity */
extern int POPULATION;
extern int WATER;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" -T  %d >> Number of times to be computed\n", TIMES);
	      printf(" -S  %f >> Step Size or time interval\n", STEP_SIZE);
	      printf(" -E  %f >> Temporal Accuracy\n", EPSILON);
	      printf(" -R  %d >> No REPLICAS or Stochastic Realizations\n", Realizations);
	      printf(" -L  %d >> Simulation Length\n", Simlength);
	      printf(" -m  %g >> Natural Mortality Rate\n", mu);
	      printf("  Recovery Rate:\n");
	      printf(" -g1 %g >> \ta1=%g\n", g1,g1);
	      printf("  Induced Mortality:\n");
	      printf(" -a1 %g >> \ta1=%g\n", a1,a1);
	      printf("  Dispersal Rate into the Water body:\n");
	      printf(" -r1 %g >> \tr1=%g\n", r1,r1);
	      printf("  Decaying Rate in Water:\n"); 
	      printf(" -o1 %g >> \to1=%g\n", o1,o1);
	      printf("  Direct Transmission Rates:\n"); 
	      printf(" -h1 %g >> \tbh1=%g\n", bh1,bh1);
	      printf("  Water-borne Transmission Rates:\n"); 
	      printf(" -w1 %g >> \tbw1=%g\n", bw1,bw1);
	      printf(" -P  %d >> Population Size in the Villages\n",POPULATION);
	      printf(" -W  %d >> Reservoir Capacity\n",WATER);
	      printf(" -f  %g >> Fraction of Water Contamination\n",f);
	      printf("\n");
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp," -T  %d >> Number of times to be computed\n", TIMES);
	      fprintf(fp," -S  %f >> Step Size or time interval\n", STEP_SIZE);
	      fprintf(fp," -E  %f >> Accuracy\n", EPSILON);
	      fprintf(fp," -R  %d >> Stochastic Realizations\n", Realizations);
	      fprintf(fp," -L  %d >> Simulation Length\n", Simlength);
	      fprintf(fp," -m  %g >> Natural Mortality Rate\n", mu);
	      fprintf(fp,"  Recovery Rate:\n");
	      fprintf(fp," -g1 %g >> \ta1=%g\n", g1,g1);
	      fprintf(fp,"  Induced Mortality:\n");
	      fprintf(fp," -a1 %g >> \ta1=%g\n", a1,a1);
	      fprintf(fp,"  Dispersal Rate into the Water body:\n");
	      fprintf(fp," -r1 %g >> \tr1=%g\n", r1,r1);
	      fprintf(fp,"  Decaying Rate in Water:\n"); 
	      fprintf(fp," -o1 %g >> \to1=%g\n", o1,o1);
	      fprintf(fp,"  Direct Transmission Rates:\n"); 
	      fprintf(fp," -h1 %g >> \tbh1=%g\n", bh1,bh1);
	      fprintf(fp,"  Water-borne Transmission Rates:\n"); 
	      fprintf(fp," -w1 %g >> \tbw1=%g\n", bw1,bw1);
	      fprintf(fp," -P  %d >> Population Size in the Villages\n",POPULATION);
	      fprintf(fp," -W  %d >> Reservoir Capacity\n",WATER);
	      fprintf(fp," -f  %g >> Fraction of Water Contamination\n",f);
  fclose(fp);
}	




