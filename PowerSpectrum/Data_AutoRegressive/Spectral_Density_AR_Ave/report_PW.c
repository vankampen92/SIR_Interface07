/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	                                                                     */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern int ovrlap;
extern int m;
extern int REPLICAS;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" Command line arguments (and default values):\n");
	      printf(" -0 %d >> Overlap option\n",ovrlap);
	      printf(" -R %d >> No. of Points in the power spectrum\n", m);
	      printf(" -m %d >> No. of REPLICAS\n", REPLICAS);
	      printf(" As an example,\n"); 
	      printf("        >> PWS -m 1000 -P 21000\n\n"); 
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp," Command line arguments (and default values):\n");
	      fprintf(fp," -0 %d >> Overlap option\n",ovrlap);
	      fprintf(fp," -R %d >> No. of REPLICAS\n", REPLICAS);
	      fprintf(fp," -m %d >> No. of Points in the power spectrum\n", m);
  fclose(fp);
}	




