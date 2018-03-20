#include "../SIR.h"
#include "../SIR_Analytic_General.h"

float reactivity(ParamSet *P)
{ 
  float **h, **v;
  float *d;
  float react;
  int nRot;

  h = matrix(1,2,1,2);
  v = matrix(1,2,1,2);
  d = vector(1,2); 

  h[1][1] = a_11(P);
  h[1][2] = (a_12(P) + a_21(P))/2.;
  h[2][1] = (a_21(P) + a_12(P))/2.;
  h[2][2] = a_22(P);

  jacobi(h, 2, d, v, &nRot);
  
  if(d[1] > d[2])
    react = d[1];
  else
    react = d[2];
  
  free_vector(d, 1,2);
  free_matrix(h, 1,2,1,2); free_matrix(v, 1,2,1,2);
  
  return (react);
}
