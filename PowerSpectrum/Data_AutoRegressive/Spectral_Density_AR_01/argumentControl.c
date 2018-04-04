#include "../../pwsp.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern char *File;
extern int ovrlap;
extern int m;

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
	    case 'O': /* Ovrlap option */
	      sscanf(argv[argcount+1],"%d", &ovrlap);
	      skip++;
	      break;
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':
	      printf(" Command line arguments (and default values):\n");
	      printf(" -h     >> help\n");
	      printf(" -0 %d >> Overlap option\n",ovrlap);
	      printf(" -m %d >> No. of Points in the power spectrum\n", m);
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




