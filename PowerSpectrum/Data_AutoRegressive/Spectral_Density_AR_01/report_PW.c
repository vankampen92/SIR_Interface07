/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	                                                                     */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern char *File;
extern int ovrlap;
extern int m;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" Command line arguments (and default values):\n");
	      printf(" -0 %d >> Overlap option\n",ovrlap);
	      printf(" -m %d >> No. of Points in the power spectrum\n", m);
	      printf(" -f  %s >> Name of the data file to be read.\n", File);
	      printf(" As an example,\n"); 
	      printf("        >> PWS -m 1000 -P 21000\n\n"); 
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp," Command line arguments (and default values):\n");
	      fprintf(fp," -0 %d >> Overlap option\n",ovrlap);
	      fprintf(fp," -m %d >> No. of Points in the power spectrum\n", m);
	      fprintf(fp," -f  %s >> Name of the data file to be read.\n", File);
	      
  fclose(fp);
}	




