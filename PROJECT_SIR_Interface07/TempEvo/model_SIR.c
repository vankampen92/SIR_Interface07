/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../SIR.h"

int SIR( Parameter_Table * Table, double Time_Initial, double * Time_Final, int * New )
{
  /* (*New) counts the number of infection events */ 
  int k;
  double Max_Probability, Total_Rate;
  double time, inter_event_time; 
  int bool_time, Event;
  
  int no_Patch               = Table->no_Patch;
  if( no_Patch > 1 ) {
    printf("This function is designed for a mono-patch system\n");
    printf("For a multi-patch system, another function should be used instead\n");
  }
  assert( no_Patch == 1 );

  Community * Village        = Table->Patch[0];   
  ParamSet * Set             = Table->Model_Parameters;
  
  int TIME_DEPENDENT_PARAM   = Table->TIME_DEPENDENCE;
  int DISCARDING_EXTINCTIONS = Table->DISCARDING_EXTINCTIONS;
  
  double time_Factor         = Table->Time_Scale_Unit;
  double ** Yp               = Table->Populations_Per_Patch;
  int *  Sum_I               = Table->Total_Populations;
  
  time = (* Time_Final);
  bool_time = 0.; (*New) = 0;
  
  while(bool_time == 0){
     
    if (TIME_DEPENDENT_PARAM > 0){ 
      evolving_Parameters(time, time_Factor, Set, TIME_DEPENDENT_PARAM);
    }
      
    Temporal_Dynamics(Village, Set, &Total_Rate, &Max_Probability, Sum_I);
 
    if(Total_Rate == 0){
      bool_time = 1;
      inter_event_time = 0.;
    }
    else{
      inter_event_time = -1./Total_Rate * log(drand48());
    }

    time += inter_event_time;
    
    if(time < (*Time_Final) && bool_time == 0){

      /* BEGIN : Stochastic Dynamic is actually performed : Village is updated accoundingly */
      Event = execute_Step(Village, Set, &Total_Rate, Max_Probability, Sum_I);
      /*   END : Stochasctic Dynamics * * * * * */

      if(Event == 0) (*New)++; /* Counting infection events... */

      if(DISCARDING_EXTINCTIONS == 1){
	if(Sum_I[1] == 0){
	  for(k=0; k < 3; k++) Yp[0][k] = (double)Sum_I[k];
	  printf("Number of infective: %d\tAn extintion has occurred\n", Sum_I[1]);
	  return(0);
	}
      }
    }
    else {
      bool_time = 1;
    }	
  }
 
  for(k=0; k<=2; k++) Yp[0][k] = Sum_I[k];
  (*Time_Final) = time - inter_event_time;

  return(0);
}
  
  





