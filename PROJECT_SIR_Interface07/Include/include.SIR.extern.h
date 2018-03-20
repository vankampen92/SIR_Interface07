/* Global Shared parameters main Program <---> ArgumentControl() */
/* Defintion Stochastic Simulation Parameters */
extern int Realizations;
extern int Simlength;
extern int No_of_Points;
extern double EPSILON;
extern int TIMES;     /* Number of times to be analyzed */
extern double STEP_SIZE; 
extern double t_0, t_1;

/* Definition Initial Conditions */
extern int I_0, M_0;      /* Inital number of infective Individuals */

/* Definition Model Parameters */
extern int No_of_Villages;
extern int POPULATION;
extern double a,d,g,mu,b,Imm;
extern double time_Factor;

/* Definition External Forcing */
extern double Per;  /* Period of the seasonal Forcing: usually 365 days */
extern double b_0, b_1;

/* Simulation Control Parameters */
extern int TIME_DEPENDENT_PARAM;
extern int DISCARDING_EXTINCTIONS;

