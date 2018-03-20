#include "../../../SIR.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern double a;
extern double mu;
extern double b, b_m, b_M;
extern double g;
extern double d, d_0, d_1;
extern double Imm;
extern int POPULATION;
extern int No_of_Points;
extern float t_f;
extern int NUMBER_of_PLOTS;
extern int SCAN_WHAT;
extern int STO_DET;
extern int TIME_DEPENDENT_PARAM;
extern int MSD;

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
	    case 'A': /* Tipus of calculation for the amplification  */
	      sscanf(argv[argcount+1],"%d", &MSD);
	      skip++;
	      break;
	    case 'S': /* Seasonal Forcing...   */
	      sscanf(argv[argcount+1],"%d", &TIME_DEPENDENT_PARAM);
	      skip++;
	      break;
	    case 'F': /* Number elements in the vector to analize 
			 different plots at different levels of the 
		         third paramter
		      */
	      sscanf(argv[argcount+1],"%d", &NUMBER_of_PLOTS);
	      skip++;
	      break;  
	    case 'W': /* Scan What? */
	      sscanf(argv[argcount+1],"%d", &SCAN_WHAT);
	      skip++;
	      break;  
	    case 'H': /* Derministic or Stochastic Dynamics */
	      sscanf(argv[argcount+1],"%d", &STO_DET);
	      skip++;
	      break;  
	    case 'N': /* Number of Points in the Transect */
	      sscanf(argv[argcount+1],"%d", &No_of_Points);
	      skip++;
	      break;
	    case 'P': /* Population Size */
	      sscanf(argv[argcount+1],"%d", &POPULATION);
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
	    case 'T': /* Backgound Environmental Transmission Rate */
	      sscanf(argv[argcount+1],"%f",&t_f);
	      skip++;
	      break;
	    case 'B': /* Initial values: Maximum and Minimum Transmission Rates  */
	      if(argv[argcount][2]=='m')
		sscanf(argv[argcount+1],"%lf",&b_m);
	      else if(argv[argcount][2]=='M')
		sscanf(argv[argcount+1],"%lf",&b_M);
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
	      printf(" -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      printf("        S = 0 (No Seasonal Forcing)\n");
	      printf("        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      printf("        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      printf("        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      printf("        where Term is +1 for school terms and -1 for holidays\n");
	      printf(" -H %d >> Stochastic (-H 0) or Derministic (-H 1)\n", STO_DET);
	      printf(" -A %d >> Stochastic Amplification:\n", MSD);
	      printf(" -A %d >> Averaging maxima (-A 0) or by the mean squared deviation(-A 1)\n",
		     MSD);
	      printf(" -W %d >> What to Scan (-W 0,  Delta Transect; -W 1, Beta Transect)\n",
		     SCAN_WHAT);
	      printf(" -N %d >> Transect scaning: %d data point\n", No_of_Points, No_of_Points);
	      printf(" -F %d >> Number of Transects\n", NUMBER_of_PLOTS);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf("\n");
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
              printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = m,M\n");
              printf("       >> b_m = %g and b_M= %g\n", b_m, b_M);
	      printf(" -D[k]   >> Maximum and Minimum Relative Infectious phase (d/g) to be scanned\n");
	      printf("       >> k = 0,1\n");
              printf("       >> d_0 = %g and d_1= %g\n", d_0, d_1);
              printf(" -mu %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -e %g >> Immgration Induced transmission rate\n", Imm);
	      printf(" -T %g >> Final time (numerical integration)\n", t_f);
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




