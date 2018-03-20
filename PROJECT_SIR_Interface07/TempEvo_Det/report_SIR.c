#include "../SIR.h"
extern float x_i, x_s;
extern float p_1, p_2;
extern double a;
extern double g;
extern double mu;
extern double b;
extern double b_0, b_1;
extern double Per;
extern double d;
extern double Imm;
extern int TIME_DEPENDENT_PARAM;
extern int POPULATION;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      printf("        S = 0 (No Seasonal Forcing)\n");
	      printf("        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      printf("        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      printf("        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      printf("        where Term is +1 for school terms and -1 for holidays\n");
	      printf(" -ti %f >> Initial Time\n", x_i);
	      printf(" -tf %f >> Final Time\n", x_s);
	      printf(" -p1 %f >> Initial fraction of susceptible\n", p_1);
	      printf(" -p2 %f >> Initial fraction of infected\n", p_2);
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" -Y %g >> Characteristic Period of the Seasonal Forcing\n",
		     Per);
	      printf("  External Transmission Rates (Immigration):\n"); 
	      printf(" -e %g\n", Imm);
	      printf(" -mu %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp, " -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      fprintf(fp, "        S = 0 (No Seasonal Forcing)\n");
	      fprintf(fp, "        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      fprintf(fp, "        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      fprintf(fp, "        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      fprintf(fp, "        where Term is +1 for school terms and -1 for holidays\n");
	      fprintf(fp," -ti %f >> Initial Time\n", x_i);
	      fprintf(fp," -tf %f >> Final Time\n", x_s);
	      fprintf(fp," -p1 %f >> Initial fraction of susceptible\n", p_1);
	      fprintf(fp," -p2 %f >> Initial fraction of infected\n", p_2);
	      fprintf(fp, " -d %g >> Natural Mortality Rate\n", d);
	      fprintf(fp, "  Induced Mortality:\n");
	      fprintf(fp, " -a %g >> \ta=%g\n", a,a);
	      fprintf(fp, "  Direct Transmission Rates (Beta):\n"); 
	      fprintf(fp, " -b %g >> \tb=%g\n", b,b);
	      fprintf(fp, " -B[k]   >> Maximum and Minimum Transmission rates\n");
	      fprintf(fp, "       >> k = 0,1\n");
	      fprintf(fp, "       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      fprintf(fp, " -Y %g >> Characteristic Period of the Seasonal Forcing\n",
		     Per);
	      fprintf(fp, "  External Transmission Rates (Immigration):\n"); 
	      fprintf(fp, " -e %g\n", Imm);
	      fprintf(fp, " -mu %g >> Loss of Immnunity Rate\n", mu);
	      fprintf(fp, 
		      " -g %g >> Recovery Rate from the infection phase\n", g);
	      fprintf(fp, " -P %d >> Population Size in the Villages\n",POPULATION);
  fclose(fp);
}	




