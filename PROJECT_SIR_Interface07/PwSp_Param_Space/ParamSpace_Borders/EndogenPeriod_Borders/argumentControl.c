#include "../../../SIR.h" 

/* Global Shared parameters main Program <---> ArgumentControl() */
extern double a;
extern double mu;
extern double b, b_0, b_1;
extern double g;
extern double d, d_0, d_1;
extern double Imm;
extern int POPULATION;
extern int No_of_Points;
extern float TOLERANCE;
extern int NUMBER_of_PLOTS;
extern int No_of_PERIOD_VALUES;
extern int SCAN_WHAT;

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
	    case 'F': /* Number of levels of immigration to plot */
	      sscanf(argv[argcount+1],"%d", &NUMBER_of_PLOTS);
	      skip++;
	      break;
	    case 'V': /* Number of Equi-Period curves */
	      sscanf(argv[argcount+1],"%d", &No_of_PERIOD_VALUES);
	      skip++;
	      break;
	    case 'H': /* What to scan */
	      sscanf(argv[argcount+1],"%d", &SCAN_WHAT);
	      skip++;
	      break;
	    case 'N': /* Number of Points in the Time Series */
	      sscanf(argv[argcount+1],"%d", &No_of_Points);
	      skip++;
	      break;
	    case 'P': /* Population Size */
	      sscanf(argv[argcount+1],"%d", &POPULATION);
	      skip++;
	      break;
	    case 'T': /* Tolerance */
	      sscanf(argv[argcount+1],"%f", &TOLERANCE);
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
	      printf("Open SIR withouth Seasonal Forcing:\n");
	      printf(" -H %d >> What to Scan (-H 0, Sucseptible; -H 1, Infective)\n",
		     SCAN_WHAT);
	      printf(" -F %d >> Number of immigration level (Different Contourn Plots)\n", 
	             NUMBER_of_PLOTS);
	      printf(" -V %d >> Number of Period values\n", No_of_PERIOD_VALUES);
	      printf(" -T %f >> Tolerance\n", TOLERANCE);
	      printf(" -N %d >> No. of Points to be saved\n", No_of_Points);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf("\n");
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
              printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
              printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" -D[k] >> Maximum and Minimum Relative Infectious phase (d/g) to be scanned\n");
	      printf("       >> k = 0,1\n");
              printf("       >> d_0 = %g and d_1= %g\n", d_0, d_1);
              printf(" -mu %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -e %g >> Immgration Induced transmission rate\n", Imm);
	      printf("\n"); 
	      printf(" As an example,\n"); 
	      printf("        >> SIR -m 0.1 -R 100000 -d 0.01 -b 20.\n\n");
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




