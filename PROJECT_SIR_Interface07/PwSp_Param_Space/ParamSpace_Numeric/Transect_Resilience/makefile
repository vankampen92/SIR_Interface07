#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*	                                                                      */
#/*                            Power Spectrum Analysis                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
HOME = /homes/dalonso
C_DIR = ${HOME}/C_pro/
ODE = ${C_DIR}/Differential_Equations/ODEqs.c 
PARAM = ${C_DIR}Stochastics/SIR/PwSp_Param_Space/ParamSpace_Numeric/
ODE_01 = ${C_DIR}/Differential_Equations/ODEqs_01.c 
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -I${C_DIR} -g #-Wall 
OBJECTS = stat.o argumentControl.o report_SIR.o main_SIR.o func_SIR.o func_Open_SIR.o func_PwSp_Analytic.o init_SIR.o resilience.o 
CC=gcc  

all: APWS

APWS:  	$(OBJECTS)   
	$(CC) -o APWS $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}stat.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}stat.c $(CXXFLAGS)

main_SIR.o: main_SIR.c ../../../SIR.h ../../../SIR_Analytic.h
	$(CC) -c main_SIR.c $(CXXFLAGS)

argumentControl.o: argumentControl.c ../../../SIR.h
	${CC} -c argumentControl.c $(CXXFLAGS)

report_SIR.o: report_SIR.c ../../../SIR.h
	${CC} -c report_SIR.c $(CXXFLAGS)

func_SIR.o:	../../../func_SIR.c ../../../SIR.h
	$(CC) -c ../../../func_SIR.c $(CXXFLAGS) 

resilience.o: ../../resilience.c ../../../SIR.h ../../../SIR_Analytic_General.h
	$(CC) -c ../../resilience.c $(CXXFLAGS) 

func_Open_SIR.o: ../../func_Open_SIR.c ../../../SIR.h ../../../SIR_Analytic_General.h
	$(CC) -c ../../func_Open_SIR.c $(CXXFLAGS) 

func_PwSp_Analytic.o: ../../func_PwSp_Analytic.c ../../../SIR.h ../../../SIR_Analytic_General.h
	$(CC) -c ../../func_PwSp_Analytic.c $(CXXFLAGS) 

init_SIR.o : ../../../init_SIR.c ../../../SIR.h
	$(CC) -c ../../../init_SIR.c $(CXXFLAGS) 

clean:  
	rm -f *.o APWS *~ core         







