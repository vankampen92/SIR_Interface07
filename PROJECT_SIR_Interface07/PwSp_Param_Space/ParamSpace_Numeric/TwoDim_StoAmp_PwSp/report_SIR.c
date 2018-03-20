#include "../../../SIR.h" 
extern double a;
extern double g;
extern double mu;
extern double b, b_0, b_1;
extern double d, d_0, d_1;
extern double Imm;
extern int No_of_Points;
extern int COHERENCE;
extern float r;
extern int NUMBER_of_PLOTS;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" -N %d >> No. of Points in the grid: %d x %d\n", 
		     No_of_Points,No_of_Points,No_of_Points); 
	      printf(" -F %d >> Number elements in the vector to analize different\n", NUMBER_of_PLOTS);
	      printf("          plots at different levels of the third paramter\n");
	      printf(" -C %d >> Control parameter for COHERENCE calculation\n", 
		     COHERENCE);
	      printf("       >> COHERENCE = 0, for calculation of overall amplification\n");
	      printf("       >> COHERENCE = 1, for calculation of relative amplification\n");
	      printf("       >> COHERENCE = 2, for calculation of both\n");
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
	      printf(" -r %f >> Coherence calculation between f_p*(1-r) and f_p*(1+r)\n", r);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp, " MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp, " -N %d >> No. of Points in the grid: %d x %d\n", 
		      No_of_Points,No_of_Points,No_of_Points);
	      fprintf(fp, " -F %d >> Number elements in the vector to analize different\n", NUMBER_of_PLOTS);
	      fprintf(fp, "          plots at different levels of the third paramter\n");
	      fprintf(fp, " -C %d >> Control parameter for COHERENCE calculation\n", 
		     COHERENCE);
	      fprintf(fp, "       >> COHERENCE = 0, for calculation of overall amplification\n");
	      fprintf(fp, "       >> COHERENCE = 1, for calculation of relative amplification\n");
	      fprintf(fp, "       >> COHERENCE = 2, for calculation of both\n");
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
	      fprintf(fp, " -r %f >> Coherence calculation between f_p*(1-r) and f_p*(1+r)\n", r);
  fclose(fp);
}	






