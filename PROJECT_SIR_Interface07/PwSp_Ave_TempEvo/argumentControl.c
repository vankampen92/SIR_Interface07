#include "../SIR.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern int No_of_Points;
extern int k_Sets;
extern int TIMES;     /* Number of times to be analyzed */
extern int TRANSIENT;
extern float STEP_SIZE; 
extern float EPSILON;
extern int Realizations;
extern int Simlength;
extern float t_0,t_1;
extern double a;
extern double mu;
extern double b;
extern double b_0;
extern double b_1;
extern double Per;  /* Period of the seasonal Forcing: usually 365 days */
extern double g;
extern double d;
extern double Imm;
extern int POPULATION;
extern int I_0, M_0;  
extern int No_of_Points;
extern int TIME_DEPENDENT_PARAM;
extern int DISCARDING_EXTINCTIONS;
extern int EQ_BETA;

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
	    case 'S': /* Seasonal Forcing...   */
	      sscanf(argv[argcount+1],"%d", &TIME_DEPENDENT_PARAM);
	      skip++;
	      break;
	    case 'D': /* Seasonal Forcing...   */
	      sscanf(argv[argcount+1],"%d", &DISCARDING_EXTINCTIONS);
	      skip++;
	      break;
	    case 'Q': /* Seasonal Forcing...   */
	      sscanf(argv[argcount+1],"%d", &EQ_BETA);
	      skip++;
	      break;
	    case 'N': /* Number of Points in the Time Series */
	      sscanf(argv[argcount+1],"%d", &No_of_Points);
	      skip++;
	      break;
	    case 'K': /* Number of segmets to be done */
	      sscanf(argv[argcount+1],"%d", &k_Sets);
	      skip++;
	      break;
	    case 'P': /* Population Size */
	      sscanf(argv[argcount+1],"%d", &POPULATION);
	      skip++;
	      break;
	    case 'R': /* Number of Realizations */
	      sscanf(argv[argcount+1],"%d", &Realizations);
	      skip++;
	      break;
	    case 'L': /* Simulation Lengnth */
	      sscanf(argv[argcount+1],"%d", &Simlength);
	      skip++;
	      break;
	    case 'I': /* Initial Number of Infective Individuals */
	      sscanf(argv[argcount+1],"%d", &I_0);
	      skip++;
	      break;
	    case 'M': /* Initial Number of Infective Individuals */
	      sscanf(argv[argcount+1],"%d", &M_0);
	      skip++;
	      break;
	    case 'T': /* Number of Times to be saved */
	      sscanf(argv[argcount+1],"%d", &TIMES);
	      skip++;
	      break;
	    case 'W': /* Number of Times to be discard */
	      sscanf(argv[argcount+1],"%d", &TRANSIENT);
	      skip++;
	      break;
	    case 't': /* Initial and final time */
	      if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%f",&t_0);
	      else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%f",&t_1);
	      else exit(0);
	      skip++;
	      break;
	    case 'Z': /* Step Size */
	      sscanf(argv[argcount+1],"%f", &STEP_SIZE);
	      skip++;
	      break;
	    case 'E': /* Error in Step Size */
	      sscanf(argv[argcount+1],"%f", &EPSILON);
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
	    case 'Y': /* Seasonal Forcing Period Per = 365 d. = 1 y. */
	      sscanf(argv[argcount+1],"%lf",&Per);
	      skip++;
	      break;
	    case 'B': /* Maximum and Minimum Transmission Rate */
	      if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%lf",&b_0);
	      else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%lf",&b_1);
	      else exit(0);
	      skip++;
	      break;
	    case 'e': /* Backgound Environmental Transmission Rate */
	      sscanf(argv[argcount+1],"%lf",&Imm);
	      skip++;
	      break;
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':
	      printf(" Command line arguments (and default values):\n");
	      printf(" -h     >> help\n");
	      printf(" -S %d  >> Seasonal Forcing Control:\n",
		     TIME_DEPENDENT_PARAM);
	      printf(" -D %d  >> Disarding time series where there is an extinction event:\n",
		     DISCARDING_EXTINCTIONS);
	      printf("      D = 0 Time Series where there is some extinctions and reintroduccions are considered\n");
	      printf("      D = 1 Only Time Series without extinctions are considered\n");
	      printf(" -Q %d  >> Seasonal Forcing Control:\n", EQ_BETA);
	      printf("      Q = 0 (Non-activation of an equivalent Transmission Rate)\n");
	      printf("      Q = 1 (Non Seasonal Forced dynamics with an equivalent Beta)\n");
	      printf("      S = 0 (No Seasonal Forcing)\n");
	      printf("      S = 1 (Sinusoidal Seasonal Forcing)\n");
	      printf("      S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      printf("      S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      printf("      where Term is +1 for school terms and -1 for holidays\n");
	      printf(" -N %d >> No. of Points to be saved\n", No_of_Points);
	      printf(" -K %d >> No. of segments in which the time series will be fragmented\n", 
		     k_Sets);
	      printf(" -R %d >> No REPLICAS or Stochastic Realizations\n", 
		     Realizations);
	      printf(" -P %d >> Population Size in the Villages\n",POPULATION);
	      printf(" -I %d >> Initial number of Infective Individuals\n",I_0);
	      printf(" -M %d >> Initial number of Immune Individuals\n",M_0);
	      printf("\n");
	      printf(" -L %d >> Simulation Length\n", Simlength);
	      printf(" -T %d >> No. of Times to be saved\n", TIMES);
	      printf(" -W %d >> No. of Times to be disarded in the transients\n", TRANSIENT);
	      printf(" -Z %f >> Step Size\n", STEP_SIZE);
	      printf(" -E %f >> Maximum Error in Step Size\n", EPSILON);
	      printf(" Times defining the Stationary Regim:\n");
	      printf(" -t0 %f >> Initial Time\n", t_0);
	      printf(" -t1 %f >> Final Time\n", t_1);
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf("  External Transmission Rates (Immigration):\n"); 
	      printf(" -e %g\n", Imm);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" -Y %g >> Characteristic Period of the Seasonal Forcing\n",
		     Per);
	      printf(" -m %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
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




