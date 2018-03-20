/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*                      Computing Deterministic Amplification                */
/*                             David Alonso, 2005 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../SIR.h"
#include "../../SIR_Analytic_General.h"

void the_Largest_Maxima(double *Maxim, int no, int no_of_Max)
{
  /* It fishes the first no_of_Max-th largest numbers from the 
     vector Maxim[] and put them at the begining of the array:
     Maxim[0] > Maxim[1] > ... > Maxim[No_of_Max-1] without
     touching the rest of points */
  int i, j, j_max;
  double Max;

  for(i=0; i<no_of_Max; i++){
    Max = Maxim[i];
    j_max = i;
    for(j=i; j<no; j++){ 
      if(Maxim[j] > Max) j_max = j;  
      Max = MAX(Max, Maxim[j]);
    }
    Maxim[j_max] = Maxim[i];
    Maxim[i] = Max;
  }
}
