/* This file, when necessary, must be include always after SIR.h */

/* These functions are defined in func_PwSp_Analytic.c */
double power_Spectrum(ParamSet *P, int I, double w);
double Alpha(ParamSet *P, int I);
double Beta(ParamSet *P, int I);
double Resonance_Frequency(ParamSet *P);
double Resonance_Frequency_Peak(ParamSet *P, int I);
double Resonance_Frequency_Peak_Approx(ParamSet *P, int I);
double Damping_Frequency(ParamSet *P);
double Maximum_Power(ParamSet *P, int I);
double Omega_2(ParamSet *P);
double Gamma_2(ParamSet *P);
double Trace(ParamSet *P);
int Condition(ParamSet *P, double factor);
int Condition_Stability(ParamSet *P);
int Exact_Condition_Peak(ParamSet *P, int I);
void Stability_General(ParamSet *P);
void eigen_Values(ParamSet *P);
double overall_Analytic(ParamSet *P, int I);
double coherence_Analytic(ParamSet *P, int I, double r, double nu_p);

/* These functions are defined in a specific file (p.e, func_Open_SIR.c)
   where the specific expressions for the values of the fix point and
   the entries of the stability matrix, A=(a_{ij}), and the correlation 
   matrix B=(b_{ij}) are given */
double Fi(ParamSet *P);
double Psi(ParamSet *P);
void Fixed_Points_General(ParamSet *P, double *S, double *I);
double a_12(ParamSet *P);
double a_21(ParamSet *P);
double a_22(ParamSet *P);
double a_11(ParamSet *P);
double b_11(ParamSet *P);
double b_12(ParamSet *P);
double b_21(ParamSet *P);
double b_22(ParamSet *P);

/* These functions are intended to calculates stability-instability borders
   and other borders in the paremeter space. They are defined in 
   func_Instability_Borders_Gen.c
*/
double Delta_2(ParamSet *P);
double Delta_4(ParamSet *P);
float zbrent_Param(float (*func)(ParamSet *P, float, float), ParamSet *P, 
		      float x1, float x2, float tol, float NS);
float border_beta_delta(float (*Function)(ParamSet *, float, float), ParamSet *P, 
			float b0, float b1, float Tolerance, float d);
float Determinant(ParamSet *P, float b, float d);
float Discriminant_2(ParamSet *P, float b, float d);
float Discriminant_4(ParamSet *P, float b, float d);

/* Functions in beta_2_4.c */
float beta_2(float d, float signe);
float beta_4(float d, float signe);

/* These functions are intended to calculate equi-period borders
   and other borders in the paremeter space. They are defined in 
   zbrent_Period.c
*/
float endogen_Period(ParamSet *P, float b, float d, float Period, int SI);
float border_Period(float (*Function)(ParamSet *, float, float, float, int), 
		    ParamSet *P, float b0, float b1, float Pe, int SI,
		    float Tolerance, float d);
float zbrent_Endogen(float (*func)(ParamSet *P, float, float, float, int), ParamSet *P,
		     float x1, float x2, float Pe, int SI, float tol, float NS);


/* Required function to calculate analytical stochastic 
   overall amplification and coherence stored in overall_Amplification.c */
void overall_Amplification_Coherence(ParamSet *P, int SI, double f_semi, double f_peak,
				     float *over, float *cohe);
double overall_Amplification(ParamSet *, int);
double coherence_value(ParamSet *, int, double, double);
double coherence_semi_analytic(ParamSet *, int, double, double);
double coherence_value_Simple(double, double);
float qromb_Accuracy(float (*func)(float), float a, float b, float EPS);
float qtrap_EPS(float (*func)(float), float a, float b, float EPS);

/* Functions to calculate the reactivity in reactivity.c */
float reactivity(ParamSet *P);

/* Functions to calculate the resilience in resilience.c */
float resilience(ParamSet *P);
float Real_Eigen_Value(ParamSet *P);
float Imaginary_Eigen_Value(ParamSet *P);

