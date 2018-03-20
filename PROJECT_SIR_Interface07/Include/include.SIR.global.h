/* Global Shared parameters main Program <---> ArgumentControl() */
/* Defintion Stochastic Simulation Parameters */
int Realizations;
int Simlength;
int No_of_Points;
double EPSILON;
int TIMES;     /* Number of times to be analyzed */
double STEP_SIZE; 
double t_0, t_1;

/* Definition Initial Conditions */
int I_0, M_0;      /* Inital number of infective Individuals */

/* Definition Model Parameters */
int No_of_Villages;
int POPULATION;
double a,d,g,mu,b,Imm;
double time_Factor;

/* Definition External Forcing */
double Per;  /* Period of the seasonal Forcing: usually 365 days */
double b_0, b_1;

/* Simulation Control Parameters */
int TIME_DEPENDENT_PARAM;
int DISCARDING_EXTINCTIONS;

