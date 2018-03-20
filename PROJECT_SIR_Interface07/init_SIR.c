#include "SIR.h"

extern double d;
extern double a;
extern double g;
extern double mu;
extern double b;
extern double Imm;
extern int POPULATION;

void Init_SumPop(Community *pVil, int *Sum_I, int Sp)
{
  int i,k;
  Community *P;
  int sum;
  
  for(k=0; k<=Sp; k++){
    sum = 0; P = pVil;
    for(i=0; i<No_of_Villages; i++, P++) sum += P->I[k];
    Sum_I[k] = sum;
  }
}

void settingParameterStruct(ParamSet *Par)
{
  Par->Beta = b;
  Par->Alpha= a;
  Par->Gamma= g;
  Par->Delta= d;
  Par->Mu = mu;
  Par->Imm = Imm;
}

void re_settingParamStruct(ParamSet *Par, double d_Rel, double b_Rel)
{
  Par->Beta = b_Rel;
  Par->Delta= d_Rel;
}

void re_setting_Equivalent_BETA(ParamSet *Par, double b_0, double b_1)
{
  double R_0, logR_0, school_days, holidays;
  
  school_days = 176.;
  holidays = 189.;

  logR_0 = log(b_0) - log(Par->Gamma) + (school_days - holidays)/365. * log(1. + b_1);

  R_0 = exp(logR_0);
  
  /* Equivalent Beta */
  Par->Beta = Par->Gamma *R_0;
}

void changingTimeScale(ParamSet *Par, double TimeScale)
{
  Par->Beta  *= TimeScale;
  Par->Alpha *= TimeScale;
  Par->Gamma *= TimeScale;
  Par->Delta *= TimeScale;
  Par->Mu    *= TimeScale;
  Par->Imm   *= TimeScale;
}
  
void Initial_Condition(Community *Village, int I_0, int M_0)
{
  Community *P;
  int i;
  
  P=Village;
  for(i=0; i<No_of_Villages; i++, P++){
    P->N = POPULATION;
    P->m = POPULATION - I_0 - M_0;
    P->n = I_0;
    P->I[0] = POPULATION - I_0 - M_0;
    P->I[1] = I_0;
    P->I[2] = M_0;
    P->position.x = 0.; P->position.y = 0.;
  }
}







