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
extern int NUMBER_of_PLOTS;
extern int SCAN_WHAT;
extern int COHERENCE;
extern double f_semiband;

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
	    case 'W': /* Scan What? */
	      sscanf(argv[argcount+1],"%d", &SCAN_WHAT);
	      skip++;
	      break; 
	    case 'C': /* Controling coherence calculation */
	      sscanf(argv[argcount+1],"%d", &COHERENCE);
	      skip++;
	      break;
	    case 's': /* Recovery from the infective phase */
	      sscanf(argv[argcount+1],"%lf", &f_semiband);
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
	      printf(" -W %d >> What to Scan (-W 0,  Delta Transect; -W 1, Beta Transect)\n",
		     SCAN_WHAT);
	      printf(" -C %d >> Control parameter for COHERENCE calculation\n", 
		     COHERENCE);
	      printf("       >> COHERENCE = 0, for calculation of overall amplification\n");
	      printf("       >> COHERENCE = 1, for calculation of relative amplification\n");
	      printf("       >> COHERENCE = 2, for calculation of both\n");
	      printf(" -s %g >> Semiband for coherence calculation\n", f_semiband);
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




