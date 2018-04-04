/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	                                                                     */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern char *File;
extern int ovrlap;
extern int m;
extern int k;
extern int No_of_Data;
extern int label;
extern float time_Step;

void modelReport(char *FileReport)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" Command line arguments (and default values):\n");
	      printf(" -0 %d >> Overlap option\n",ovrlap);
	      printf(" -D %d >> No. of Data Points in the file to be read\n", No_of_Data);
	      printf(" -m %d >> No. of Points in the power spectrum\n", m);
	      printf(" -k %d >> No. of Segments to divide the time series\n", k);
	      printf(" -T %f  >> Sampling time interval\n",time_Step);
	      printf(" -l %d  >> File label\n", label);
	      printf(" -f  %s >> Name of the data file to be read.\n", File);
	      printf(" As an example,\n"); 
	      printf("        >> PWS -m 1000 -P 21000\n\n"); 
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(FileReport, "w");
	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp," Command line arguments (and default values):\n");
	      fprintf(fp," -D %d >> No. of Data Points in the file to be read\n", No_of_Data);
	      fprintf(fp," -0 %d >> Overlap option\n",ovrlap);
	      fprintf(fp," -m %d >> No. of Points in the power spectrum\n", m);
	      fprintf(fp," -k %d >> No. of Segments to divide the time series\n", k);
	      fprintf(fp," -T %f  >> Sampling time interval\n",time_Step);
	      fprintf(fp," -l %d  >> File label\n", label);
	      fprintf(fp," -f  %s >> Name of the data file to be read.\n", File);
	      
  fclose(fp);
}	




