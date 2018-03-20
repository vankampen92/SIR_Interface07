#include "../../../SIR.h" 
extern double a;
extern double g;
extern double mu;
extern double b, b_m, b_M;
extern double d, d_0, d_1;
extern double Imm;
extern int POPULATION;
extern int No_of_Points;
extern int NUMBER_of_PLOTS;
extern int SCAN_WHAT;
extern int COHERENCE;
extern double f_semiband;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" -W %d >> What to Scan (-W 0,  Delta Transect; -W 1, Beta Transect)\n",
		     SCAN_WHAT);
	      printf(" -F %d >> Number of Transects\n", NUMBER_of_PLOTS);
	      printf(" -N %d >> Transect scaning: %d data point\n", No_of_Points, No_of_Points);
	      printf(" -C %d >> Control parameter for COHERENCE calculation\n", 
		     COHERENCE);
	      printf("       >> COHERENCE = 0, for calculation of overall amplification\n");
	      printf("       >> COHERENCE = 1, for calculation of relative amplification\n");
	      printf("       >> COHERENCE = 2, for calculation of both\n");
	      printf(" -s %g >> Semiband for coherence calculation\n", f_semiband);
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
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp, " MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp, " -W %d >> What to Scan (-W 0,  Delta Transect; -W 1, Beta Transect)\n",
		     SCAN_WHAT);
	      fprintf(fp, " -F %d >> Number of Transects\n", NUMBER_of_PLOTS);
	      fprintf(fp, " -N %d >> Transect scaning: %d data point\n", 
		      No_of_Points, No_of_Points);
	      fprintf(fp, " -C %d >> Control parameter for COHERENCE calculation\n", 
		     COHERENCE);
	      fprintf(fp, "       >> COHERENCE = 0, for calculation of overall amplification\n");
	      fprintf(fp, "       >> COHERENCE = 1, for calculation of relative amplification\n");
	      fprintf(fp, "       >> COHERENCE = 2, for calculation of both\n");
	      fprintf(fp, " -s %g >> Semiband for coherence calculation\n", f_semiband);
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
  fclose(fp);
}	





