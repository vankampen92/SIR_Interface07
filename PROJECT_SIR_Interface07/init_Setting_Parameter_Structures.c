#include "./Include/SIR.h"

#include "./Include/include.SIR.extern.h"

void setting_Dynamic_Model_Parameter_Structure( ParamSet * P)
{
  P->Beta = b;
  P->Alpha= a;
  P->Gamma= g;
  P->Delta= d;
  P->Mu = mu;
  P->Imm = Imm;
}
  
void setting_Parameter_Table(ParamSet * P, Community ** PATCH, Parameter_Table * Table, 
			     int no_Types)
{
  int i;

  Table->Model_Parameters = P;
  Table->Patch            = PATCH;

  Table->POPULATION       = POPULATION;
  
  Table->TIME_DEPENDENCE       = TIME_DEPENDENT_PARAM;
  Table->DICARDING_EXTINCTIONS = DISCARDING_EXTINCTIONS;
  
  Table->no_Patch              = No_of_Villages;

  Table->Simulation_Length     = Simlength;
  
  Table->Time_0                = t_0;         
  Table->Infectious_Time_0     = I_0;     
  Table->Susceptible_Time_0    = M_0;   
  Table->Time_Scale_Unit       = time_Factor;

  Table->Time                  = t_0;

  Table->Populations_Per_Patch = (double**)calloc(No_of_Villages, sizeof(double *) );
  for(i=0; i<No_of_Villages; i++) {
    Table->Populations_Per_Patch[i] = (double *)calloc(no_Types, sizeof(double) );
  }

  Table->Total_Populations = (int *)calloc( no_Types, sizeof(int) );
}

void free_Parameter_Table( Parameter_Table * Table )
{
  int i;

  for(i=0; i<No_of_Villages; i++) {
    free (Table->Populations_Per_Patch[i]);
  }
  free (Table->Populations_Per_Patch);
  
  free (Table->Total_Populations);
}

void  Initial_Condition_System( Parameter_Table * Table)
{
  Community * Village;
  int i;
  int n, I_0, M_0; 

  n   = Table->no_Patch;

  N   = Table->POPULATION;
  I_0 = Table->Infectious_Time_0;
  M_0 = Table->Susceptible_Time_0;
  
  for(i; i < n; i++){
    Initial_Condition_Individual_Patch( Table->PATCH[i], S_0, I_0, N );
  }  

  Table->Populations_Per_Patch[0] = (double)M_0;
  Table->Populations_Per_Patch[1] = (double)I_0;
  Table->Populations_Per_Patch[2] = (double)N - (double)M_0 - (double)I_0;

  if( n == 1 ){
    Table->Total_Populations[0] = M_0;
    Table->Total_Populations[1] = I_0;
    Table->Total_Populations[2] = N - M_0 - I_0;
  }
  else{
    printf(" Attention!!! This code only works for one globally mixed community.\n");
    printf(" The program will exit...");
    /* Init_SumPop(Village, Sum_I, 2); */
    exit(0);
  }
}





