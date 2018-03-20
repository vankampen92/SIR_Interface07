#include "../../../SIR.h" 
extern double a;
extern double g;
extern double mu;
extern double b, b_m, b_M;
extern double d, d_0, d_1;
extern double Imm;
extern int POPULATION;
extern int No_of_Points;
extern float t_f;
extern int NUMBER_of_PLOTS;
extern int SCAN_WHAT;
extern int STO_DET;
extern int TIME_DEPENDENT_PARAM;
extern int MSD;

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
	      printf(" -H %d >> Stochastic (-H 0) or Derministic (-H 1)\n", STO_DET);
	      printf(" -A %d >> Stochastic Amplification:\n", MSD);
	      printf(" -A %d >> Averaging maxima (-A 0) or by the mean squared deviation(-A 1)\n",
		     MSD);
	      printf(" -W %d >> What to Scan (-W 0,  Delta Transect; -W 1, Beta Transect)\n",
		     SCAN_WHAT);
	      printf(" -F %d >> Number of Transects\n", NUMBER_of_PLOTS);
	      printf(" -N %d >> Transect scaning: %d data point\n", No_of_Points, No_of_Points);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
              printf("       >> b_m = %g and b_M= %g\n", b_m, b_M);
	      printf(" -D[k]   >> Maximum and Minimum Relative Infectious phase (d/g) to be scanned\n");
	      printf("       >> k = 0,1\n");
              printf("       >> d_0 = %g and d_1= %g\n", d_0, d_1);
	      printf(" -m %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -e %g >> Immgration Induced transmission rate\n", Imm);
	      printf(" -T %f >> Immgration Induced transmission rate\n", t_f);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp, " MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp, " -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      fprintf(fp, "        S = 0 (No Seasonal Forcing)\n");
	      fprintf(fp, "        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      fprintf(fp, "        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      fprintf(fp, "        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      fprintf(fp, "        where Term is +1 for school terms and -1 for holidays\n");
	      fprintf(fp, " -H %d >> Stochastic (-H 0) or Derministic (-H 1)\n", STO_DET);
	      fprintf(fp, " -A %d >> Stochastic Amplification:\n", MSD);
	      fprintf(fp, " -A %d >> Averaging maxima (-A 0) or by the mean squared deviation(-A 1)\n",
		     MSD);
	      fprintf(fp, " -W %d >> What to Scan (-W 0,  Delta Transect; -W 1, Beta Transect)\n",
		     SCAN_WHAT);
	      fprintf(fp, " -F %d >> Number of Transects\n", NUMBER_of_PLOTS);
	      fprintf(fp, " -N %d >> Transect scaning: %d data point\n", 
		      No_of_Points, No_of_Points);
	      fprintf(fp, " -P %d >> Population Size in the Villages\n",POPULATION);
	      fprintf(fp, " -d %g >> Natural Mortality Rate\n", d);
	      fprintf(fp, "  Induced Mortality:\n");
	      fprintf(fp, " -a %g >> \ta=%g\n", a,a);
	      fprintf(fp, "  Direct Transmission Rates (Beta):\n"); 
	      fprintf(fp, " -b %g >> \tb=%g\n", b,b);
	      fprintf(fp, " -B[k]   >> Maximum and Minimum Transmission rates\n");
	      fprintf(fp, "       >> k = 0,1\n");
              fprintf(fp, "       >> b_m = %g and b_M= %g\n", b_m, b_M);
	      fprintf(fp, " -D[k]   >> Maximum and Minimum Relative Infectious phase (d/g) to be scanned\n");
	      fprintf(fp, "       >> k = 0,1\n");
              fprintf(fp, "       >> d_0 = %g and d_1= %g\n", d_0, d_1);
	      fprintf(fp, " -m %g >> Loss of Immnunity Rate\n", mu);
	      fprintf(fp, " -g %g >> Recovery Rate from the infection phase\n", g);
	      fprintf(fp, " -e %g >> Immgration Induced transmission rate\n", Imm);
	      fprintf(fp, " -T %f >> Immgration Induced transmission rate\n", t_f);
  fclose(fp);
}	






