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

extern int TIME_DEPENDENT_PARAM;
extern int M_0; 
extern float time_Factor;

void SIR_StocRea(Community *Village, ParamSet *Set, 
		 int Simlength, int modul, double t_final, 
		 int Sum_I[], double *Time, int **Yp, int I_0)
{
  int i,ii,iii, k,Event;
  double Max_Probability, Total_Rate, time, inter_event_time; 
  static int replicaCounter = 0;
  char Prefix[No_of_Villages][14];
  
  Initial_Condition(Village, I_0, M_0);
  Init_SumPop(Village, Sum_I, 2);
  
  time = 0.;
  inter_event_time = 1.; /* it must be positive, so that Temporal dynamics() 
			    may iniciate the first total rate computation!! */
  iii  = 0; 
  for (i=0;i<Simlength;i++){

    if(inter_event_time > 0.){
      if (TIME_DEPENDENT_PARAM > 0){ 
	evolving_Parameters(time, time_Factor, Set, TIME_DEPENDENT_PARAM);
      }
      Temporal_Dynamics(Village, Set, &Total_Rate, &Max_Probability, Sum_I);
    }
    inter_event_time = -1./Total_Rate * log(drand48());   
    time += inter_event_time; 
    Event = execute_Step(Village, Set, &Total_Rate, Max_Probability, Sum_I);
    ii=i%modul;
    if(ii==0){
      Time[iii] = time;
      for(k=0; k<=2; k++) Yp[iii][k] = Sum_I[k];
      iii++;
      ShowVillage_StrainEvolution(Village, 0, time);
      //SaveVillageContent("pFr", Village, 1, replicaCounter, time);
      //if(time > 5.) /*Saving Stationary Time Series */
      //SaveTimeSeries(Village,replicaCounter, time);
      //Press_Key();
    }
  }
  /* End of simulation loop */      
  replicaCounter++;
  printf("New Replica... No = %d\n", replicaCounter);
}
  
  





