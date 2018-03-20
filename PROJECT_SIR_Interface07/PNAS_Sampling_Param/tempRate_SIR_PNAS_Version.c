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

void Temporal_Dynamics(Community *Village, ParamSet *pa, 
		       double *Rate, double *max_Probability, int *I_sum)
{
  int i,k;
  Community *P;
  
  P = Village;
  *max_Probability = 0.; *Rate = 0.;
  /* Transitions allowed within the patch system */
  for(i=0; i<No_of_Villages; i++, P++){
  
    P->rate[0] = 0.;
    /* Infection rate of one susceptible through human contact or 
       environmental immigration */
    P->rate[0]+=pa->Beta*P->I[1]/(double)P->N + pa->Imm;
    
    P->H_infect = P->rate[0];
    
    /* Recovery plus Death-Birth (Constant Population) rate of Infected 
       individuals */
    P->rate[1] = pa->Delta + pa->Alpha + pa->Gamma; 
    
    /* Death-Birth (Constant population) of Recovered Individuals */
    P->rate[2] = pa->Mu + pa->Delta;

    /* Vilage Transition Rate */  
    P->ratePatch = 0.;
    for(k=0; k<=2; k++) {
      P->ratePatch += (double)P->I[k] * P->rate[k];
      P->rToI[k] = (double)P->I[k] * P->rate[k];
    }
    *Rate += P->ratePatch;
    (*max_Probability) = MAX(*max_Probability, P->ratePatch);
    
    if(*Rate < 0.){
      printf("\n");
      printf("R, Total rate for an eventual change in system configuration\n");
      printf("R = %g\n", *Rate);
      printf("As R is zero or negative, no change is possible\n");
      printf("\n");
      exit(0);
    }
  }
}







