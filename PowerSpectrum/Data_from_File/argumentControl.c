#include "../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern char *File;
extern int ovrlap;
extern int No_of_Data;
extern int m;
extern int k;
extern int label;
extern float time_Step;

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
	    case 'f': /* Name of the file to be read */
	      sscanf(argv[argcount+1],"%s", File);
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
	    case 'D': /* Number of Data */
	      sscanf(argv[argcount+1],"%d", &No_of_Data);
	      skip++;
	      break;
	    case 'T': /* Number of Data */
	      sscanf(argv[argcount+1],"%f", &time_Step);
	      skip++;
	      break;
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':
	      printf(" Command line arguments (and default values):\n");
	      printf(" -h     >> help\n");
	      printf(" -0 %d  >> Overlap option\n",ovrlap);
	      printf(" -D %d  >> No. of Data Points in the file to be read\n", No_of_Data);
	      printf(" -m %d  >> No. of Points in the power spectrum\n", m);
	      printf(" -k %d  >> No. of Segments to divide the time series\n", k);
	      printf(" -T %f  >> Sampling time interval\n",time_Step);
	      printf(" -l %d  >> File label\n", label);
	      printf(" -f  %s >> Name of the data file to be read.\n", File);
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




