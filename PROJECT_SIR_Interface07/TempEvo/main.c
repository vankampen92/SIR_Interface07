/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../Include/SIR.h"

#include "../Include/include.SIR.global.h"

void fill_FileNames(char Pref[2][25], int);

void fill_FileNames(char Pref[2][25], int No)
{
  int i;
  for(i=0; i < No; i++){
    Pref[i][0]='\0';   
    fitxer(Pref[i], "sp_", i, ".dat");
  }
}

int main(int argc, char **argv)
{
  /* Temporal evolution of the number of strains in each stochastic realization */
  int i,j,k,modul, spNumber, Bad_Times; 
  /* Time_values[TIMES] = 
     {0.2, 0.3, 0.4, 0.5, 0.7, 1., 2., 5., 6., 8., 10., 12., 15., 20.}; */

  double ave, var;
  double Psi,Fi;
 
  char Files[2][25];
  FILE *fp[2]; 
  FILE *Fp; 
  char File[12];

  ParamSet Dynamic_Model_Par;  
  Parameter_Table Model_Par_Table;
  
#include "../Include/include.SIR.default.c"  
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);
  
  /* BEGIN : Memmory allocation * * * * * * * * * * * * * * * * */
  Communty ** PATCH = (Community **)malloc( No_of_Villages * sizeof( Village *) );
  for( i = 0; i<No_of_Villages; i++ ){    PATCH[i] = Patch__Memmory__Allocation(3);   }

  double * Time_values = (double *)calloc( TIMES, sizeof(double) );
  
  double ** summ      = (double **)calloc( TIMES, sizeof(double *) );
  for (i=0; i<TIMES; i++)  summ[i]   = (double *)calloc( 3, sizeof(double) );

  double ** summ_var   = (double **)calloc( TIMES, sizeof(double *) );
  for (i=0; i<TIMES; i++)  summ_var[i]   = (double *)calloc( 3, sizeof(double) );

  int * count = (int *)calloc(TIMES, sizeof(int) );
    
  int ** New_Infections   = (int **)calloc( Realizations, sizeof(double *) );
  for (i=0; i < Realizations; i++)  New_Infections[i]   = (int *)calloc( TIMES, sizeof(int) );
  /*   END : Memmory allocation * * * * * * * * * * * * * * * * */

  /* BEGIN : Initial Setting * * * * * * * * * * * * * * * */
  setting_Dynamic_Model_Parameter_Structure(&Dynamic_Model_Par);
  
  setting_Parameter_Table( &Dynamic_Model_Par, &PATCH, &Model_Par_Table, 3 );

  Fixed_Points_General(&Dynamic_Model_Par, &Psi, &Fi);
  I_0 = POPULATION*Fi; M_0 = POPULATION*(1.-Fi-Psi);
  /* Alternatively: */
  //I_0 = 0.0001 * POPULATION; M_0 = POPULATION - 0.1 * POPULATION - I_0;
  
  for(i=0; i<TIMES; i++) Time_values[i] = t_0 + (double)(i+1) * STEP_SIZE;
  printf("Sampling Interval (days): %g\n", STEP_SIZE);
  /*  END : Initial Setting * * * * * * * * * * * * * * * * */

  printf("Entering Generation of Stochastic Realizations...\n");
  Press_Key();
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */  
  /* Main loop: a number of REALIZATIONS (stochastic temporal evolutions) is computed */
  for (i=0;i<Realizations;i++){
    /* Input variables: */
    /* 
       . Model_Par_Table, a comprehensive model parameter table (see definition in SIR.h) 
       . Time_values[], sampling times, a vector of 'TIMES' different times */
       
    /* Output variables: */
    /* . New Infections is the number of accumulated infections per time interval.
         This is the usual measured variable in a hospital */
    /* . summ and (summ_var) are the sum (and the summ of squares) of each variable 
         type (S, I, R) at each time across simulations.
         They will allow to calculate averages and variances */
    /* . Bad_Times is a measure of the performance of the sampling frequency. 
         If Bad_Times is high, interval times should be choosen smaller */
 
    S_T_O_C_H_A_S_T_I_C___T_E_M_P_O_R_A_L___E_V_O_L_U_T_I_O_N(i, 
							      &Model_Par_Table, 
							      TIMES, Time_values,
							      summ, summ_var, count,
							      New_Infections,
							      Bad_Times);
    
    /* End of the i-th STOCHASTIC REALIZATIONS */
    printf("Realization: %d of a total of %d\n", i+1, Realizations);
    printf("Time failed in %d occasions in %d time steps\n", Bad_Times, TIMES);
    printf("EPSILON might be too small!\n");
    printf("Try to choose a larger EPSILON -E %g\n",EPSILON);
  }
  /* End of STOCHASTIC REALIZATIONS */
  /* Saving Realizations */
  
  /* Re-scaling times (annual)*/
  for(i=0; i<TIMES; i++)
    Time_values[i] = t_0 + (double)(i+1) * STEP_SIZE * 1./Per;
  
  for(i=0; i < Realizations; i++){
    File[0]='\0';  fitxer(File, "new", i, ".dat"); Fp = fopen(File, "w");
    for(j=0; j < TIMES; j++){  
      fprintf(Fp,"%g\t%d\n", Time_values[j], New_Infections[i][j]);
    }
    fclose(Fp);
  }
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  
  /* BEGIN : Saving Results * * * * * * * * * * * * * * * * */
  fill_FileNames(Files, 2);
  /*Opening files */
  for(i=0; i<=1; i++) fp[i]=fopen(Files[i],"w");
  /* Re-scaling times (dayly)*/
  for(i=0; i<TIMES; i++)
    Time_values[i] = t_0 + (double)(i+1) * STEP_SIZE;

  for (j=0; j<TIMES; j++){
    /* Saving the strain temporal evolution */
    for(k=0; k<=1; k++){
      if(count[j] > 0){
	ave = summ[j][k]/(float)count[j];
	var = summ_var[j][k]/(float)count[j] - ave*ave;
	if(var >= 0)
	  var = sqrt(var);
	else
	  var = 0.;
	fprintf(fp[k],"%g\t%g\t%g\n",Time_values[j],ave,var);
      }
    }
  }
  /*Closing files */
  for(i=0; i<=1; i++) fclose(fp[i]);
 
  modelReport("report.txt");
  /* END : Saving Results * * * * * * * * * * * * * * * * * */

  /* BEGIN : Freeing All Memmory * * * * * * * * * * * * * * */
  for( i = 0; i<No_of_Villages; i++ ){ Patch__Memmory__Free( PATCH[i] ); }
  free (PATCH);

  free (Time_values);

  for( i=0; i<TIMES; i++) {
    free( summ[i] );
    free( summ_var[i] );
  }
  free(summ); 
  free(summ_var);

  free (count);
  
  for(i=0; i<Realizations; i++) {
    free( New_Infections[i] );
  }
  free( New_Infections );

  free_Parameter_Table( Model_Par_Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */
 
  printf("\nEnd of progam\n");
  return (0);
}





