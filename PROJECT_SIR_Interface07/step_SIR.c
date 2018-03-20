/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "SIR.h"

int execute_Step(Community *Pop, ParamSet *pa, 
		 double *Total_Rate, double max_Probability,
		 int Sum_I[]) 
{
  /* 
     Here a rejection method is implemented to choose the event to occur.
     This involves that the Village within the system is elected at random, and
     that the Total_Probability is in fact the largest transition rate for any
     species within the whole patch system.
  */
  int x, I;
  double P1;
  Community *pVil;
  int Event;

  /* Hierarchic procedure to find the even to occur... */  
  /* The event occurs within the villages... */
  if(No_of_Villages == 1) x = 0;
  else                    x = which_Village(max_Probability, Pop);
  pVil = &Pop[x];  
  I = Discrete_Sampling(pVil->rToI, 3) - 1; 
  /* I = 0, 1, 2 */
  if(I == 0){                      
    /* The event is affecting the Susceptible population: Infection Event */
    Sum_I[1]++; pVil->I[1]++; pVil->n++; 
    Sum_I[0]--; pVil->I[0]--; pVil->m--; 
  }
  else if(I==1){
    P1 = pa->Gamma/pVil->rate[1]; 
    if(drand48() < P1){
      /* The event is affecting the Infected population: Recovery Event */
      Sum_I[1]--;  pVil->I[1]--; pVil->n--; 
      Sum_I[2]++;  pVil->I[2]++;  
    }
    else{
      /* The event is affecting the Infected population: Birth-Death Event */
      Sum_I[1]--; pVil->I[1]--; pVil->n--; 
      Sum_I[0]++; pVil->I[0]++; pVil->m++; 
    }
  }
  else if(I==2){
    /* The event is affecting on already recovered individual: Death-ReBirth */
    Sum_I[2]--; pVil->I[2]--; 
    Sum_I[0]++; pVil->I[0]++; pVil->m++; 
  }
  else{
    /* Something wrong!!! */
    printf("I should be between either 0, 1 or 2\n");
    printf("I = %d\n", I);
    Press_Key();
  }
  Event = I;
  return (Event);
}






