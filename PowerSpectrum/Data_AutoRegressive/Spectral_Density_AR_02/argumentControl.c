#include "../../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern int No_of_Points;
extern int REPLICAS;
extern int ONE;

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
	    case 'P': /* Number of Points in the Time Series */
	      sscanf(argv[argcount+1],"%d", &No_of_Points);
	      skip++;
	      break;
	    case 'R': /* Number of Replicas */
	      sscanf(argv[argcount+1],"%d", &REPLICAS);
	      skip++;
	      break;
	    case 'O': /* Number of Replicas */
	      sscanf(argv[argcount+1],"%d", &ONE);
	      skip++;
	      break;
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':
	      printf(" Command line arguments (and default values):\n");
	      printf(" -h     >> help\n");
	      
	      printf(" -P %d >> No. of Points in the time series\n",No_of_Points);
	      printf(" -0 %d >> Contral parameter:\n",ONE);
	      printf("       >> ONE = 1, dfour1(...) is used\n");
	      printf("       >> ONE = 0, twofft(...) is used\n");
	      printf(" -R %d >> No. of Averaged Replicas in the power spectrum\n", REPLICAS);
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




