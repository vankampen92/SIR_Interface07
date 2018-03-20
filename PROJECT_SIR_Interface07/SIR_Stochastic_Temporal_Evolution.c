#include "./Include/SIR.h"

void S_T_O_C_H_A_S_T_I_C___T_E_M_P_O_R_A_L___E_V_O_L_U_T_I_O_N(int i, 
							       Parameter_Table * Table, 
							       int TIMES, double * Time_values, double EPSILON,
							       double * summ, double * summ_var, 
							       int * count,
							       int ** New_Infections,
							       int * Bad_Times)
{
  FILE *FP; char file[12];
  int j, k;
  int New; /* Accumulate number of infections */
  double Time_Initial, Time_Final;
  int N;

  N = Table->POPULATION;
  
  /* Each stochastic realization will be saved in a different file */
  file[0]='\0';  fitxer(file, "re_", i, ".dat"); FP = fopen(file, "w");
  
  /* BEGIN : Initial Conditions */
  Initial_Condition_System( Table );
  /* will have to call Init_SumPop(Village, Sum_I, 2); */
  /*   END : Initial Conditions */
  
  (*Bad_Times) = 0;
  for (j=0; j < TIMES-1; j++){
    
    Time_Initial = Time_values[j];
    Time_Final = Time_values[j+1];
    
    /* Now the system will be advance from Time_Initial to Time_Final */
    SIR( Table, Time_Initial, &Time_Final, &new );
    
    if(Time_Final > Time_values[j+1] - EPSILON){
      //printf("New Infections: %d\n", new);
      //printf("Actual Infectious: %d\t%d\n\n", Strain[1], Village[0].I[1]);
      
      for(k=0; k < 3; k++){
	summ[j][k] += Table->Total_Populations[k]/(float)N;
	summ_var[j][k] += Table->Total_Populations[k]/(float)N * Table->Total_Populations[k]/(float)N;	
      }
      count[j]++;
    }
    else{
      (*Bad_Times)++;
      printf("%d\t%g\t%g\n",j, Time_Final, Time_values[j]);
    }
    
    fprintf(FP,"%g\t%d\t%d\n", Time_Final, Table->Total_Populations[0], Table->Total_Populations[1]);
    New_Infections[i][j] = new;
    /* Time[0] is always the last time which is the closest to  Time_values[j] */
    /* go further to the next time */
  }
  fclose(FP);
}
