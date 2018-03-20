#include "../../SIR.h" 
extern int Realizations;
extern int Simlength;
extern double a;
extern double g;
extern double mu;
extern double Imm;
extern double b;
extern double b_0;
extern double b_1;
extern double Per;  /* Period of the seasonal Forcing: usually 365 days */
extern double d;
extern int POPULATION;
extern int TIME_DEPENDENT_PARAM;
extern float l_0, l_1;
extern float d_0, d_1;
extern float r_0, r_1;

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
	      printf(" -R %d >> No REPLICAS or Stochastic Realizations\n", 
		     Realizations);
	      printf(" -L %d >> Simulation Length\n", Simlength);
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" Ranges to explore:\n");
	      printf(" -l[k]   >> Range in Betas (Transmission Rates)\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> l_0 = %g and l_1= %g\n", l_0, l_1);
	      printf(" -D[k]   >> Range in Deltas (Tunrover rates)\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> d_0 = %g and d_1= %g\n", d_0, d_1);
	      printf(" -r[k]   >> Range in Immgration levels\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> r_0 = %g and r_1= %g\n", r_0, r_1);
	      printf(" -Y %g >> Characteristic Period of the Seasonal Forcing\n",
		     Per);
	      printf(" -mu %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf("  External Transmission Rate (Immigration):\n"); 
	      printf(" -e %g\n", Imm);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
              fprintf(fp, " -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      fprintf(fp, "        S = 0 (No Seasonal Forcing)\n");
	      fprintf(fp, "        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      fprintf(fp, "        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      fprintf(fp, "        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      fprintf(fp, "        where Term is +1 for school terms and -1 for holidays\n");
      	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp, " -R %d >> No REPLICAS or Stochastic Realizations\n", 
		      Realizations);
	      fprintf(fp, " -L %d >> Simulation Length\n", Simlength);
	      fprintf(fp, " -d %g >> Natural Mortality Rate\n", d);
	      fprintf(fp, "  Induced Mortality:\n");
	      fprintf(fp, " -a %g >> \ta=%g\n", a,a);
	      fprintf(fp, "  Direct Transmission Rates (Beta):\n"); 
	      fprintf(fp, " -b %g >> \tb=%g\n", b,b);
              fprintf(fp, " -B[k]   >> Maximum and Minimum Transmission rates\n");
	      fprintf(fp, "       >> k = 0,1\n");
	      fprintf(fp, "       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      fprintf(fp, " Ranges to explore:\n");
	      fprintf(fp, " -l[k]   >> Range in Betas (Transmission Rates)\n");
	      fprintf(fp, "       >> k = 0,1\n");
	      fprintf(fp, "       >> l_0 = %g and l_1= %g\n", l_0, l_1);
	      fprintf(fp, " -D[k]   >> Range in Deltas (Tunrover rates)\n");
	      fprintf(fp, "       >> k = 0,1\n");
	      fprintf(fp, "       >> d_0 = %g and d_1= %g\n", d_0, d_1);
	      fprintf(fp, " -r[k]   >> Range in Immgration levels\n");
	      fprintf(fp, "       >> k = 0,1\n");
	      fprintf(fp, "       >> r_0 = %g and r_1= %g\n", r_0, r_1);	      
	      fprintf(fp, " -Y %g >> Characteristic Period of Seasonal Forcing\n",
		     Per);  
	      fprintf(fp, " -mu %g >> Loss of Immnunity Rate\n", mu);
	      fprintf(fp, 
		      " -g %g >> Recovery Rate from the infection phase\n", g);
	      fprintf(fp, "  External Transmission Rate (Immigration):\n"); 
	      fprintf(fp, " -e %g\n", Imm);
              fprintf(fp, 
		      " -P %d >> Population Size in the Villages\n",POPULATION);
  fclose(fp);
}	




