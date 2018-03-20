#include "./Include/SIR.h"

Community * Patch__Memmory__Allocation( int no_Types )
{
  Community * P = (Community *)malloc(sizeof ( Community ));
  
  P->N = 0;             /* Total Patch capacity: Human Constant Population     */
  P->m = 0;             /* Total human population of susceptibles m = I[0]     */
  P->n = 0;             /* n = I[1]  Number of infected                        */
  P->I = (int *)calloc( no_Types, sizeof(int) );         
  
  P->rate = (double *)calloc( no_Types, sizeof(double) ); 
  P->rToI = (double *)calloc( no_Types, sizeof(double) ); 
  
  P->ratePatch = 0.0;  /* Transition probability at this patch                    */
  P->H_infect  = 0.0;
  
  //struct point position;  /* Patch Geographical Coordinates  */ 

  return(P);
}

void Initial_Condition_Individual_Patch( Community * P, int S_0, int I_0, int N )
{
  P->N    = N;
  P->I[0] = S_0;
  P->I[1] = I_0;
  P->I[2] = N - S_0 - I_0;
}

void Patch__Memmory__Free( Community * P )
{
  free( P->rate );
  free( P->rToI );
  free( P->I    );
}

