/* My GSL definitions for common utilities */
#ifndef THE_BASIC_C_HEADER_FILES
#define THE_BASIC_C_HEADER_FILES
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#endif

#ifndef GSL_HEADER_FILES
#define GSL_HEADER_FILES
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#endif

void GSL_Init_Random_Seed(const gsl_rng *);
void GSL_Init_Random_Seed_from_File(const gsl_rng *);

int Discrete_Sampling_Rejection_Method(double p_max, double * a, int NoEvents);
int Discrete_Sampling_Old(double * a, int NoEvents);
int Discrete_Sampling(double *a, int NoEvents);
int Discrete_Sampling_Cummulative(double rate, double *R_A_T_E, int NoEvents);
int Discret_Sampling_High_Performance(double rate, double *R_A_T_E, int NoEvents); 
  
int da_gsl_vector_isnan(const gsl_vector *);
int da_vector_isnan( double * , int );

int fprintf_to_File_Matrix_gsl(FILE * , gsl_matrix * , int, int); 
