#include "../../../SIR.h" 
extern double a;
extern double g;
extern double mu;
extern double b, b_0, b_1;
extern double d, d_0, d_1;
extern double Imm;
extern int POPULATION;
extern int No_of_Points;
extern float TOLERANCE;
extern int NUMBER_of_PLOTS;
extern int No_of_PERIOD_VALUES;
extern int SCAN_WHAT;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf("Open SIR withouth Seasonal Forcing:\n");
	      printf(" -H %d >> What to Scan (-H 0, Sucseptible; -H 1, Infective)\n",
		     SCAN_WHAT);
	      printf(" -F %d >> Number of Contourns\n", NUMBER_of_PLOTS);
	      printf(" -V %d >> Number of Period values\n", No_of_PERIOD_VALUES);
	      printf(" -T %f >> Tolerance\n", TOLERANCE);
	      printf(" -N %d >> No. of Points to be saved\n", No_of_Points);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf("where f_r is the spectrum peak (resonance frequency.)\n");
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
              printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" -D[k]   >> Maximum and Minimum Relative Infectious phase (d/g) to be scanned\n");
	      printf("       >> k = 0,1\n");
              printf("       >> d_0 = %g and d_1= %g\n", d_0, d_1);
	      printf(" -m %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -e %g >> Immgration Induced transmission rate\n", Imm);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp, " MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp, "Open SIR withouth Seasonal Forcing:\n");
	      fprintf(fp, " -H %d >> What to Scan (-H 0, Sucseptible; -H 1, Infective)\n",
		     SCAN_WHAT);
	      fprintf(fp, " -F %d >> Number of Contourns\n", NUMBER_of_PLOTS);
	      fprintf(fp, " -V %d >> Number of Period values\n", No_of_PERIOD_VALUES);
	      fprintf(fp, " -T %f >> Tolerance\n", TOLERANCE);
	      fprintf(fp, " -N %d >> No. of Points to be saved\n", No_of_Points);
	      fprintf(fp, " -P %d >> Population Size in the Villages\n",POPULATION);
	      fprintf(fp, " -d %g >> Natural Mortality Rate\n", d);
	      fprintf(fp, "  Induced Mortality:\n");
	      fprintf(fp, " -a %g >> \ta=%g\n", a,a);
	      fprintf(fp, "  Direct Transmission Rates (Beta):\n"); 
	      fprintf(fp, " -b %g >> \tb=%g\n", b,b);
	      fprintf(fp, " -B[k]   >> Maximum and Minimum Transmission rates\n");
	      fprintf(fp, "       >> k = 0,1\n");
              fprintf(fp, "       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      fprintf(fp, " -D[k]   >> Maximum and Minimum Relative Infectious phase (d/g) to be scanned\n");
	      fprintf(fp, "       >> k = 0,1\n");
              fprintf(fp, "       >> d_0 = %g and d_1= %g\n", d_0, d_1);
	      fprintf(fp, " -m %g >> Loss of Immnunity Rate\n", mu);
	      fprintf(fp, " -g %g >> Recovery Rate from the infection phase\n", g);
	      fprintf(fp, " -e %g >> Immgration Induced transmission rate\n", Imm);
  fclose(fp);
}	






