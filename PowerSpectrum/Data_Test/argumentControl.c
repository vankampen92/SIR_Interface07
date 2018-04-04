#include "../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern double wr;
extern double noise;
extern double amplitud;
extern int ovrlap;
extern int m;
extern int k;
extern int label;
extern float time_Step;
extern int REPLICAS;

void ArgumentControl(int argc, char **argv)
{
  int argcount, skip;
  
  for(argcount=1; argcount<argc; argcount+=skip)
    {
      if(argv[argcount][0] == '-')
	{
	  skip = 1;
	  switch(argv[argcount][1])
	    {
	    case 'R': /* Number of Replicas */
	      sscanf(argv[argcount+1],"%d", &REPLICAS);
	      skip++;
	      break;
	    case 'm': /* Number of Points in the estimted Power Spectrum */
	      sscanf(argv[argcount+1],"%d", &m);
	      skip++;
	      break;
	    case 'l': /* Number of Points in the estimted Power Spectrum */
	      sscanf(argv[argcount+1],"%d", &label);
	      skip++;
	      break;
	    case 'k': /* Number of Segments */
	      sscanf(argv[argcount+1],"%d", &k);
	      skip++;
	      break;
	    case 'O': /* Ovrlap option */
	      sscanf(argv[argcount+1],"%d", &ovrlap);
	      skip++;
	      break;
	    case 'T': /* Time Step */
	      sscanf(argv[argcount+1],"%f", &time_Step);
	      skip++;
	      break;
	    case 'A': /* Amplitude */
	      sscanf(argv[argcount+1],"%lf", &amplitud);
	      skip++;
	      break;
	    case 'w': /* Endogenous Frequency */
	      sscanf(argv[argcount+1],"%lf", &wr);
	      skip++;
	      break;
	    case 'n': /* Noise */
	      sscanf(argv[argcount+1],"%lf", &noise);
	      skip++;
	      break;
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':
	      printf(" Command line arguments (and default values):\n");
	      printf(" -h     >> help\n");
	      printf(" -R %d  >> Number of Stochastic Replicas\n", REPLICAS);
	      printf(" -A %g  >> Amplitude\n",amplitud);
	      printf(" -w %g  >> Endogenous Frequency\n",wr);
	      printf(" -n %g  >> Noise\n",noise);
	      printf(" -0 %d  >> Overlap option\n",ovrlap);
	      printf(" -m %d  >> No. of Points in the power spectrum\n", m);
	      printf(" -k %d  >> No. of Segments to divide the time series\n", k);
	      printf(" -T %f  >> Sampling time interval\n",time_Step);
	      printf(" -l %d  >> File label\n", label);
	      printf(" As an example,\n"); 
	      printf("        >> PWS -m 1000 -P 21000\n\n"); 
	      exit(0);
	    }
	}
      else
	{
	  printf("**invalid command line flag >%c<\n",argv[argcount][0]);
	  printf("try -h for help.\n");
	  exit(0);
	}
    }
}	




