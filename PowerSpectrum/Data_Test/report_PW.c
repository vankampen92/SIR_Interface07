/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*	                                                                     */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern int ovrlap;
extern int m;
extern int k;
extern int label;
extern float time_Step;
extern double amplitud;
extern double wr;
extern double noise;
extern int REPLICAS;

void modelReport(char *FileReport)
{
  FILE *fp;
	      printf(" MODEL RUNNING PARAMETERS:\n");
	      printf(" -R %d  >> Number of Stochastic Replicas\n", REPLICAS);
	      printf(" -A %g  >> Amplitude\n",amplitud);
	      printf(" -w %g  >> Endogenous Frequency\n",wr);
	      printf(" -n %g  >> Noise\n",noise);
	      printf(" Command line arguments (and default values):\n");
	      printf(" -0 %d >> Overlap option\n",ovrlap);
	      printf(" -m %d >> No. of Points in the power spectrum\n", m);
	      printf(" -k %d >> No. of Segments to divide the time series\n", k);
	      printf(" -T %f  >> Sampling time interval\n",time_Step);
	      printf(" -l %d  >> File label\n", label);
	      printf(" As an example,\n"); 
	      printf("        >> PWS -m 1000 -P 21000\n\n"); 
	      printf("\n"); 
	      /* * * * * */
  fp = fopen(FileReport, "w");
	      fprintf(fp," MODEL RUNNING PARAMETERS:\n");
	      fprintf(fp," Command line arguments (and default values):\n");
	      fprintf(fp," -R %d  >> Number of Stochastic Replicas\n", REPLICAS);
	      fprintf(fp," -A %g  >> Amplitude\n",amplitud);
	      fprintf(fp," -w %g  >> Endogenous Frequency\n",wr);
	      fprintf(fp," -n %g  >> Noise\n",noise);
	      fprintf(fp," -0 %d >> Overlap option\n",ovrlap);
	      fprintf(fp," -m %d >> No. of Points in the power spectrum\n", m);
	      fprintf(fp," -k %d >> No. of Segments to divide the time series\n", k);
	      fprintf(fp," -T %f  >> Sampling time interval\n",time_Step);
	      fprintf(fp," -l %d  >> File label\n", label);	      	      
  fclose(fp);
}	




