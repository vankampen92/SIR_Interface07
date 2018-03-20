#include "../../SIR.h" 
extern double a;
extern double g;
extern double mu;
extern double b;
extern double d;
extern double Imm;
extern int POPULATION;
extern int No_of_Points;
extern double factor;
extern float r;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" -N %d >> No. of Points to be saved\n", No_of_Points);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf(" -F %g >> Calculated Spectrum from 0 to %g*f_r\n", factor,
		     factor);
	      printf("where f_r is the spectrum peak (resonance frequency.)\n");
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -m %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -e %g >> Immgration Induced transmission rate\n", Imm);
	      printf(" Coherence calculation (if it applies)\n");
	      printf(" -r %f >> Coherence calculation between f_p*(1-r) and f_p*(1+r)\n", r);
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp, " MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp, " -N %d >> No. of Points to be saved\n", No_of_Points);
	      fprintf(fp, " -P %d >> Population Size in the Villages\n",POPULATION);
	      fprintf(fp, " -d %g >> Natural Mortality Rate\n", d);
	      fprintf(fp, "  Induced Mortality:\n");
	      fprintf(fp, " -a %g >> \ta=%g\n", a,a);
	      fprintf(fp, "  Direct Transmission Rates (Beta):\n"); 
	      fprintf(fp, " -b %g >> \tb=%g\n", b,b);
	      fprintf(fp, " -m %g >> Loss of Immnunity Rate\n", mu);
	      fprintf(fp, " -g %g >> Recovery Rate from the infection phase\n", g);
	      fprintf(fp, " -e %g >> Immgration Induced transmission rate\n", Imm);
	      fprintf(fp, " Coherence calculation (if it applies!)\n");
	      fprintf(fp, " -r %f >> Coherence calculation between f_p*(1-r) and f_p*(1+r)\n", r);
  fclose(fp);
}	






