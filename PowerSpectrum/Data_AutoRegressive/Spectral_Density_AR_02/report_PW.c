/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	                                                                     */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern int No_of_Points;
extern int REPLICAS;
extern int ONE;

void modelReport(char *File)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" Command line arguments (and default values):\n");
	      printf(" -P %d >> No. of Points in the time series\n",No_of_Points);
	      printf(" -m %d >> No. of averaged REPLICAS in the power spectrum\n", REPLICAS);
	      printf(" -f  %s >> Name of the data file to be read.\n", File);
	      printf(" -P %d >> Contral parameter:\n",ONE);
	      printf("       >> ONE = 1, dfour1(...) is used\n");
	      printf("       >> ONE = 0, twofft(...) is used\n");
	      printf(" As an example,\n"); 
	      printf("        >> PWS -m 1000 -P 21000\n\n"); 
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(File, "w");
	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp," Command line arguments (and default values):\n");
	      fprintf(fp," -P %d >> No. of Points in data file\n",No_of_Points);
	      fprintf(fp," -m %d >> No. of Points in the power spectrum\n", REPLICAS);
	      fprintf(fp," -P %d >> Contral parameter:\n",ONE);
	      fprintf(fp,"       >> ONE = 1, dfour1(...) is used\n");
	      fprintf(fp,"       >> ONE = 0, twofft(...) is used\n");
  fclose(fp);
}	




