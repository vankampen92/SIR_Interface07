#include "../SIR.h" 
extern int TIMES;     /* Number of times to be analyzed */
extern float t_0,t_1;
extern float STEP_SIZE; 
extern float EPSILON;
extern int Realizations;
extern int I_0, M_0;
extern double a;
extern double g;
extern double mu;
extern double b;
extern double b_0;
extern double b_1;
extern double Imm;
extern double Per;  /* Period of the seasonal Forcing: usually 365 days */
extern double d;
extern int POPULATION;
extern int TIME_DEPENDENT_PARAM;
extern int DISCARDING_EXTINCTIONS;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" -D %d  >> Disarding time series where there is an extinction:\n",
		     DISCARDING_EXTINCTIONS);
	      printf(" -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      printf("        S = 0 (No Seasonal Forcing)\n");
	      printf("        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      printf("        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      printf("        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      printf("        where Term is +1 for school terms and -1 for holidays\n");
	      printf(" -R %d >> No REPLICAS or Stochastic Realizations\n", 
		     Realizations);
	      printf(" -T %d >> No. of Times to be saved\n", TIMES);
	      printf(" -t0 %f >> Initial Time\n", t_0);
	      printf(" -t1 %f >> Final Time\n", t_1);
	      printf(" -Z %f >> Step Size\n", STEP_SIZE);
	      printf(" -E %f >> Maximum Error in Step Size\n", EPSILON);
	      printf(" -I %d >> Initial number of Infective Individuals\n",I_0);
	      printf(" -M %d >> Initial number of Immune Individuals\n",M_0);
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" -mu %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf("  External Transmission Rate (Immigration):\n"); 
	      printf(" -e %g\n", Imm);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
           fprintf(fp, " MODEL RUNNING PARAMETERS:\n");
	   fprintf(fp," -D %d  >> Disarding time series where there is an extinction:\n",
		     DISCARDING_EXTINCTIONS);
	   fprintf(fp, " -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	   fprintf(fp, "        S = 0 (No Seasonal Forcing)\n");
	   fprintf(fp, "        S = 1 (Sinusoidal Seasonal Forcing)\n");
	   fprintf(fp, "        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	   fprintf(fp, "        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	   fprintf(fp, "        where Term is +1 for school terms and -1 for holidays\n");
	   fprintf(fp, " -R %d >> No REPLICAS or Stochastic Realizations\n", 
		  Realizations);
	   fprintf(fp, " -T %d >> No. of Times to be saved\n", TIMES);
	   fprintf(fp, " -t0 %f >> Initial Time\n", t_0);
	   fprintf(fp, " -t1 %f >> Final Time\n", t_1);
	   fprintf(fp, " -Z %f >> Step Size\n", STEP_SIZE);
	   fprintf(fp, " -E %f >> Maximum Error in Step Size\n", EPSILON);
	   fprintf(fp, " -I %d >> Initial number of Infective Individuals\n",I_0);
	   fprintf(fp, " -M %d >> Initial number of Immune Individuals\n",M_0);
	   fprintf(fp, " -d %g >> Natural Mortality Rate\n", d);
	   fprintf(fp, "  Induced Mortality:\n");
	   fprintf(fp, " -a %g >> \ta=%g\n", a,a);
	   fprintf(fp, "  Direct Transmission Rates (Beta):\n"); 
	   fprintf(fp, " -b %g >> \tb=%g\n", b,b);
	   fprintf(fp, " -B[k]   >> Maximum and Minimum Transmission rates\n");
	   fprintf(fp, "       >> k = 0,1\n");
	   fprintf(fp, "       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	   fprintf(fp, " -mu %g >> Loss of Immnunity Rate\n", mu);
	   fprintf(fp, " -g %g >> Recovery Rate from the infection phase\n", g);
	   fprintf(fp, "  External Transmission Rate (Immigration):\n"); 
	   fprintf(fp, " -e %g\n", Imm);
	   fprintf(fp, " -P %d >> Population Size in the Villages\n",POPULATION);
	   
  fclose(fp);
}	




