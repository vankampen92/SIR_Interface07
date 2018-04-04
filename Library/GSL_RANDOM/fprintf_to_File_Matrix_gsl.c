#include "GSL_stat.h"

int fprintf_to_File_Matrix_gsl(FILE * Out, gsl_matrix * A, int MM, int NN)
{
  int i,j;
  for (i=0;i<MM;i++){
    for (j=0;j<NN;j++){
      fprintf(Out,"%g\t",gsl_matrix_get(A,i,j));
    }
    fprintf(Out,"\n");
  }
  fprintf(Out,"\n");
  return 0;
}
