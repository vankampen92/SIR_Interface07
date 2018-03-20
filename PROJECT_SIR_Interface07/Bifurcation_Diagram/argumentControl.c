#include "../SIR.h" 
/* Global Shared parameters main Program <---> ArgumentControl() */
extern int No_of_Points; 
extern float x_i, x_s;
extern float p_1, p_2;
extern double a;
extern double mu;
extern double b;
extern double b_0;
extern double b_1;
extern double b_m;
extern double b_M;
extern double Per;  /* Period of the seasonal Forcing: usually 365 days */
extern double si;
extern double g;
extern double d;
extern double Imm;
extern int TIME_DEPENDENT_PARAM; 
extern int POPULATION;

extern int MODEL;
extern int No_of_Beta;
extern int No_of_Lattice;
extern int SYSTEMATIC;

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
	    case 'M': /* Random Sampling of Initial Conditions */
	      sscanf(argv[argcount+1],"%d", &MODEL);
	      skip++;
	      break;    
	    case 'R': /* Random Sampling of Initial Conditions */
	      sscanf(argv[argcount+1],"%d", &SYSTEMATIC);
	      skip++;
	      break;  
	    case 'L': /* No_of_Lattice * No_of_Lattice different initial conditions */
	      sscanf(argv[argcount+1],"%d", &No_of_Lattice);
	      skip++;
	      break;
	    case 'K': /* No of Points in the discretization of Beta */
	      sscanf(argv[argcount+1],"%d", &No_of_Beta);
	      skip++;
	      break;
	    case 'S': /* Seasonal Forcing...   */
	      sscanf(argv[argcount+1],"%d", &TIME_DEPENDENT_PARAM);
	      skip++;
	      break;
	    case 'P': /* Population in Numbers */
	      sscanf(argv[argcount+1],"%d", &POPULATION);
	      skip++;
	      break;  
	    case 'N': /* Number of Points Numerical Integration */
	      sscanf(argv[argcount+1],"%d", &No_of_Points);
	      skip++;
	      break;  
	    case 'p': /* Initial Conditions */
	      if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%f",&p_1);
	      else if(argv[argcount][2]=='2')
		sscanf(argv[argcount+1],"%f",&p_2);
	      else exit(0);
	      skip++;
	      break;
	    case 't': /* Initial and final time */
	      if(argv[argcount][2]=='i')
		sscanf(argv[argcount+1],"%f",&x_i);
	      else if(argv[argcount][2]=='f')
		sscanf(argv[argcount+1],"%f",&x_s);
	      else exit(0);
	      skip++;
	      break;
	    case 'g': /* Recovery from the infective phase */
	      sscanf(argv[argcount+1],"%lf", &g);
	      skip++;
	      break;
	    case 'x': /* Rate at which individuals leave the exposed class */
	      sscanf(argv[argcount+1],"%lf", &si);
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
	      else if(argv[argcount][2]=='m')
		sscanf(argv[argcount+1],"%lf",&b_m);
	      else if(argv[argcount][2]=='M')
		sscanf(argv[argcount+1],"%lf",&b_M);
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
	      printf(" -L %d  >> %d * %d different initial conditions\n",
		     No_of_Lattice, No_of_Lattice, No_of_Lattice);
	      printf(" -K %d  >> No of Points in the discretization of Beta\n",
		     No_of_Beta);
	      printf(" -R %d  >> Systematic (-R 0) or Random (-R 1) scaning of inital conditions\n",
		     SYSTEMATIC);
	      printf(" -M %d  >> Model to choose:\n", MODEL);
	      printf("        M = 0 SIR (constant size)\n"); 
	      printf("        M = 1 SIR (Non constant size)\n");
	      printf("        M = 2 SEIR (Non constanct size)\n");
	      printf(" -S %d  >> Seasonal Forcing Control:\n", TIME_DEPENDENT_PARAM);
	      printf("        S = 0 (No Seasonal Forcing)\n");
	      printf("        S = 1 (Sinusoidal Seasonal Forcing)\n");
	      printf("        S = 2 (beta = b_0 * (1 + b_1)^(Term)\n");
	      printf("        S = 3 (beta = b_0 * (1. + Term * b_1)\n");
	      printf("        where Term is +1 for school terms and -1 for holidays\n");
	      printf(" -P %d >> Population Size in the Villages\n", POPULATION,POPULATION);
	      printf(" -N %d >> Number of Points to be saved: %d\n", No_of_Points, No_of_Points);
	      printf(" -ti %f >> Initial Time\n", x_i);
	      printf(" -tf %f >> Final Time\n", x_s);
	      printf(" -p1 %f >> Initial fraction of susceptible\n", p_1);
	      printf(" -p2 %f >> Initial fraction of infective\n", p_2);
	      printf(" -d %g >> Natural Mortality Rate\n", d);
	      printf("  Induced Mortality:\n");
	      printf(" -a %g >> \ta=%g\n", a,a);
	      printf("  Direct Transmission Rates (Beta):\n"); 
	      printf(" -b %g >> \tb=%g\n", b,b);
	      printf(" -B[k]   >> Maximum and Minimum Transmission rates\n");
	      printf("       >> k = 0,1\n");
	      printf("       >> b_0 = %g and b_1= %g\n", b_0, b_1);
	      printf(" -B[k]   >> Maximum and Minimum values in the bifurcation diagram\n");
	      printf("       >> k = n,M\n");
	      printf("       >> b_m = %g and b_M= %g\n", b_m, b_M);
	      printf(" -Y %g >> Characteristic Period of the Seasonal Forcing\n",
		     Per);
	      printf("  External Transmission Rates (Immigration):\n"); 
	      printf(" -e %g\n", Imm);
	      printf(" -mu %g >> Loss of Immnunity Rate\n", mu);
	      printf(" -g %g >> Recovery Rate from the infection phase\n", g);
	      printf(" -x %g >> Rate at which individuals leave the exposed class\n", si);
	      printf("\n"); 
	      printf(" As an example,\n"); 
	      printf("        >> SIR -m 2.-d 0.01 -b 20.\n\n");
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



