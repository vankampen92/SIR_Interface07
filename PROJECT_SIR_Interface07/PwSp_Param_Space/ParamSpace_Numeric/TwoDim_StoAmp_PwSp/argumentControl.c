#include "../../../SIR.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern double a;
extern double mu;
extern double b, b_0, b_1;
extern double g;
extern double d, d_0, d_1;
extern double Imm;
extern int No_of_Points;
extern int COHERENCE;
extern float r;
extern int NUMBER_of_PLOTS;

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
	    case 'F': /* Number elements in the vector to analize 
			 different plots at different levels of the 
		         third paramter
		      */
	      sscanf(argv[argcount+1],"%d", &NUMBER_of_PLOTS);
	      skip++;
	      break;  
	    case 'N': /* Number of Points in the grid */
	      sscanf(argv[argcount+1],"%d", &No_of_Points);
	      skip++;
	      break;
	    case 'C': /* Number of Points in the grid */
	      sscanf(argv[argcount+1],"%d", &COHERENCE);
	      skip++;
	      break;
	    case 'r': /* Coherence calculation between f_p*(1-r) and f_p*(1+r) */
	      sscanf(argv[argcount+1],"%f", &r);
	      skip++;
	      break;
	    case 'g': /* Recovery from the infective phase */
	      sscanf(argv[argcount+1],"%lf", &g);
	      skip++;
	      break;
            case 'm': /* Loss of immunity rate */
	      sscanf(argv[argcount+1],"%lf", &mu);
	      skip++;
	      break;
	    case 'd': /* Natural Mortality */
	      sscanf(argv[argcount+1],"%lf", &d);
	      skip++;
	      break;
	    case 'a': /* Pathogen Induced Mortality */
	      sscanf(argv[argcount+1],"%lf",&a);
	      skip++;
	      break;
	    case 'b': /* Human to Human Transmission Rate  */
	      sscanf(argv[argcount+1],"%lf",&b);
	      skip++;
	      break;
	    case 'e': /* Backgound Environmental Transmission Rate */
	      sscanf(argv[argcount+1],"%lf",&Imm);
	      skip++;
	      break;
	    case 'B': /* Initial values: Maximum and Minimum Transmission Rates  */
	      if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%lf",&b_0);
	      else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%lf",&b_1);
	      else exit(0);
	      skip++;
	      break;
            case 'D': /* Initial values: Maximum and Minimum Transmission Rates  */
	      if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%lf",&d_0);
	      else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%lf",&d_1);
	      else exit(0);
	      skip++;
	      break;
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':
	      printf(" Command line arguments (and default values):\n");
	      printf(" -h    >> help\n");
	      printf(" -N %d >> No. of Points in the grid: %d x %d\n", 
		     No_of_Points,No_of_Points,No_of_Points);
	      printf(" -F %d >> Number elements in the vector to analize different\n", NUMBER_of_PLOTS);
	      printf("          plots at different levels of the third paramter\n"); 
	      printf("\n");
	      printf(" -C %d >> Control parameter for COHERENCE calculation\n", 
		     COHERENCE);
	      printf("       >> COHERENCE = 0, for calculation of overall amplification\n");
	      printf("       >> COHERENCE = 1, for calculation of relative amplification\n");
	      printf("       >> COHERENCE = 2, for calculation of both\n");	     
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
              printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
              printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" -D[k]   >> Maximum and Minimum Relative Infectious phase (d/g) to be scanned\n");
	      printf("       >> k = 0,1\n");
              printf("       >> d_0 = %g and d_1= %g\n", d_0, d_1);
              printf(" -m %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -e %g >> Immgration Induced transmission rate\n", Imm);
	      printf(" -r %f >> Coherence calculation between f_p*(1-r) and f_p*(1+r)\n", r);
	      printf("\n"); 
	      printf(" As an example,\n"); 
	      printf("        >> SIR -d 0.01 -b 20.\n\n");
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




