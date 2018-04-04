#include "../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern int No_of_Points;
extern char *File;
extern float t_0;

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
	    case 'P': /* Number of Points in the Time Series */
	      sscanf(argv[argcount+1],"%d", &No_of_Points);
	      skip++;
	      break;
	    case 'T': /* Number of Points in the Time Series */
	      sscanf(argv[argcount+1],"%f", &t_0);
	      skip++;
	      break;
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':
	      printf(" Command line arguments (and default values):\n");
	      printf(" -h     >> help\n");
	      printf(" -P %d >> No. of Points in the time series\n",No_of_Points);
	      printf(" -T %d >> Initial time after transients\n",t_0);
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




