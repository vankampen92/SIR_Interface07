#include "../SIR.h" 
extern int Realizations;
extern int Simlength;
extern double a;
extern double g;
extern double mu;
extern double b;
extern double b_0;
extern double b_1;
extern double Per;
extern double d;    
extern double Imm;
extern int POPULATION;
extern int I_0, M_0;  
extern float t_0,t_1;
extern int nP;
extern int TIME_DEPENDENT_PARAM;
extern int EQ_BETA;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL COMMAND LINE RUNNING PARAMETERS:\n");
	      printf(" -Q %d  >> Seasonal Forcing Control:\n", EQ_BETA);
	      printf("      Q = 0 (Non-activation of an equivalent Transmission Rate)\n");
	      printf("      Q = 1 (Non Seasonal Forced dynamics with an equivalent Beta)\n");
	      printf(" -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      printf("        S = 0 (No Seasonal Forcing)\n");
	      printf("        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      printf("        S = 2 (beta = b_0 * (1 + b_1)^(Term))\n");
	      printf("        S = 3 (beta = b_0 * (1. + Term * b_1))\n");
	      printf("        where Term is +1 for school terms and -1 for holidays\n");
	      printf(" -R %d >> No REPLICAS or Stochastic Realizations\n", 
		     Realizations);
	      printf(" -L %d >> Simulation Length\n", Simlength);
	      printf(" -F %d >> No. of Points entering Spectral Analysis\n", nP);
	      printf(" Times defining the Staonary Regim:\n");
	      printf(" -t0 %f >> Initial Time\n", t_0);
	      printf(" -t1 %f >> Final Time\n", t_1);
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
	      printf(" -Y %g >> Characteristic Period of the Seasonal Forcing\n",
		     Per);
	      printf("  External Transmission Rates (Immigration):\n"); 
	      printf(" -e %g\n", Imm);
	      printf(" -m %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp," MODEL COMMAND LINE RUNNING PARAMETERS:\n");
	      fprintf(fp, " -Q %d  >> Seasonal Forcing Control:\n", EQ_BETA);
	      fprintf(fp, "  Q = 0 (Non-activation of an equivalent Transmission Rate)\n");
	      fprintf(fp, "  Q = 1 (Non Seasonal Forced dynamics with an equivalent Beta)\n");
	      fprintf(fp, " -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      fprintf(fp, "        S = 0 (No Seasonal Forcing)\n");
	      fprintf(fp, "        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      fprintf(fp, "        S = 2 (beta = b_0 * (1 + b_1)^(Term))\n");
	      fprintf(fp, "        S = 3 (beta = b_0 * (1. + Term * b_1))\n");
	      fprintf(fp, "        where Term is +1 for school terms and -1 for holidays\n");
	      fprintf(fp, " -R %d >> No REPLICAS or Stochastic Realizations\n", 
		     Realizations);
	      fprintf(fp, " -L %d >> Simulation Length\n", Simlength);
	      fprintf(fp, " -F %d >> No. of Points entering Spectral Analysis\n", nP);
	      fprintf(fp, " Times defining the Staonary Regim:\n");
	      fprintf(fp, " -t0 %f >> Initial Time\n", t_0);
	      fprintf(fp, " -t1 %f >> Final Time\n", t_1);
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
	      fprintf(fp, " -Y %g >> Characteristic Period of the Seasonal Forcing\n",
		     Per);
	      fprintf(fp, "  External Transmission Rates (Immigration):\n"); 
	      fprintf(fp, " -e %g\n", Imm);
	      fprintf(fp, " -m %g >> Loss of Immnunity Rate\n", mu);
	      fprintf(fp, " -g %g >> Recovery Rate from the infection phase\n", g);
	      fprintf(fp, " -P %d >> Population Size in the Villages\n",POPULATION);
	      fprintf(fp, "\n"); 
  fclose(fp);
}	





