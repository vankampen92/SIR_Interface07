typedef struct parameterSEIRinfo
{
  double b;   /* Basal Transmision rate */
  double b_0; /* Basal Transmision rate */
  double b_1; /* Seasonality Parameter  */
  double si;  /* 1./si, average time in the Exposed, but no infected phase */
  double g ;  /* 1./g, average time in the Infectious phase */
  double nu;  /* Per capita birth rate */
  double d;   /* Per capita mortality rate */
  double a;   /* Disease-induced death rate */
  double mu;  /* Rate of loss of immunity */
  double Imm; /* Immigration parameter */
  double Per; /* Period of seasonal forcing */
}ParamSEIR;

typedef struct parameterSIRinfo
{
  double b;    /* Basal Transmision rate */
  double b_0;  /* Basal Transmision rate */
  double b_1;  /* Seasonality Parameter  */
  double g;    /* 1./g, average time in the Infectious phase */
  double nu;   /* Per capita birth rate */
  double d;   /* Per capita mortality rate */
  double a;   /* Disease-induced death rate */
  double mu;  /* Rate of loss of immunity */
  double Imm; /* Immigration parameter */
  double Per; /* Period of seasonal forcing */
}ParamSIR;

void default_values_Measles_02_SEIR(ParamSEIR *P);
void default_values_Measles_02_SIR(ParamSIR *P);
