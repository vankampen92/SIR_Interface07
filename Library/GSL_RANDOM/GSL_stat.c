#include "GSL_stat.h"

int da_vector_isnan( double *y, int n)
{
  int i, NOT_a_NUMBER;

  NOT_a_NUMBER = 0;
  
  for(i = 0; i < n; i++) {
    if ( gsl_isnan( y[i]) == 1) NOT_a_NUMBER++;
  }

  return(NOT_a_NUMBER);
}

int da_gsl_vector_isnan(const gsl_vector * y)
{
  double x;
  int i, NOT_a_NUMBER;

  NOT_a_NUMBER = 0;
  
  for(i = 0; i < y->size; i++) {
    x = gsl_vector_get(y, i);
    if ( gsl_isnan(x) == 1) NOT_a_NUMBER++;
  }

  return(NOT_a_NUMBER);
}

