#include "SIR.h" 

extern double Per;  /* Period of the seasonal Forcing: usually 365 days */
extern double b_0;
extern double b_1;

void evolving_Parameters(double time, float time_Factor, 
			 ParamSet *p, int Typus)
{
  /* time (days) = Time (dimless units) /  time_Factor; */
  /* time_Factor turns time into a time measured in days. 
     If time is already measured in days, then time_Factor should
     be 1.
     Notice that the forcing period $Per$ should be measured in the 
     same time units as times.
     Here, time_Factor is defined as with T^(-1) units.
  */
  double Term, lnBeta;
  
  switch (Typus)
    {
    case 1: /* Sinusoidal Seasonal Forcing */
      p->Beta = b_0 * (1. + b_1 * sin(2.* M_PI * time/Per));
      break;
    case 2: /* School term Seasonal Forcing (Keeling-Rohani-Grenfell, 2001) */
      Term = Time_of_the_Year(time/time_Factor);
      //Term = Time_of_the_Year_Two_Terms(time/time_Factor);
      lnBeta = log(b_0) + Term * log(1. + b_1);
      p->Beta = exp(lnBeta);   
      //printf("School Term: %g\tTransmission Rate, Beta(t=%g) = %g\n", Term, time, p->Beta);
      break;
    case 3: /* School term Seasonal Forcing (Earn-Rohani-Bolker-Grenfell, 2000) */
      Term = Time_of_the_Year(time/time_Factor);
      //Term = Time_of_the_Year_Two_Terms(time/time_Factor);
      p->Beta = b_0 * (1.+ b_1 * Term);
      //printf("School Term: %g\tTransmission Rate, Beta(t=%g) = %g\n", Term, time, p->Beta); 
      break;
    default:
      printf("Invalid Seasonal Forcing (1, 2, 3)\n");
      printf(" TIME_DEPENDENT_PARAM = %d\n", Typus);
      exit(0);
    }	
}

double Time_of_the_Year(double time)
{
  /* School Terms in England */
  /* Time should be given as an input in days */
  double rest;
  double School_Term;
  
  rest = (int)time%365;  /* The current day of the current year */
  
  if(rest >=1. && rest <= 6.)
    School_Term = -1.;  
  else if(rest >=100. && rest <= 115.)
    School_Term = -1.;
  else if(rest >=200. && rest <= 251.)
    School_Term = -1.;
  else if(rest >=300. && rest <= 307.)
    School_Term = -1.;
  else if(rest >=356. && rest <= 365.)
    School_Term = -1.;
  else
    /* During School Terms */
    School_Term = +1.;
  
  return(School_Term);
}

double Time_of_the_Year_Two_Terms(double time)
{
  /* Two Terms within a year in England */
  /* Time should be given as an input in days */
  double rest;
  double School_Term;
  
  rest = (int)time%365;  /* The current day of the current year */
  
  if(rest >=1. && rest <= 182.56)
    /* During the first half of the year */
    School_Term = -1.;  
  else
    /* During the second half of the year */
    School_Term = +1.;
  
  return(School_Term);
}
