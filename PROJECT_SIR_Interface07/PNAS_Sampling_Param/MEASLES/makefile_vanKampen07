#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/*                                  SIR MODEL                                */
#/*	                                                                      */
#/*                          (CONSTANT COMMUNITY SIZE)                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
C_DIR = /homes/dalonso/C_pro/
#MAIN = main_SIR_S0S3
MAIN = main_SIR_C_of_V
MEASLES =_Measles
ODE = ${C_DIR}/Differential_Equations/ODEqs.c 
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -I${C_DIR} #-g #-Wall -g  
OBJECTS = IO_SIR.o stat.o argumentControl.o report_SIR.o model_SIR.o func_SIR.o ${MAIN}.o init_SIR.o step_SIR.o tempRate_SIR_PNAS_Version.o evolving_Parameters.o sampling_Parameters${MEASLES}.o peak_To_trough_Ratio.o relativeAmplification_Forcing.o ODE_definition.o ODEqs.o func_PwSp_Analytic.o func_Open_SIR.o 
CC=gcc  

all: SIR

SIR:  $(OBJECTS)   
	$(CC) -o SIR $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}/stat.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}stat.c $(CXXFLAGS)

${MAIN}.o: ${MAIN}.c ../../SIR.h
	$(CC) -c ${MAIN}.c $(CXXFLAGS) 

argumentControl.o: argumentControl.c ../../SIR.h
	${CC} -c argumentControl.c $(CXXFLAGS)

model_SIR.o:	../model_SIR.c ../../SIR.h
	$(CC) -c ../model_SIR.c $(CXXFLAGS)

report_SIR.o : report_SIR.c ../../SIR.h
	$(CC) -c report_SIR.c $(CXXFLAGS) 

peak_To_trough_Ratio.o : ../../peak_To_trough_Ratio.c ../../SIR.h
	$(CC) -c ../../peak_To_trough_Ratio.c $(CXXFLAGS) 

relativeAmplification_Forcing.o : ../../relativeAmplification_Forcing.c ../../SIR.h
	$(CC) -c ../../relativeAmplification_Forcing.c $(CXXFLAGS)

sampling_Parameters${MEASLES}.o : sampling_Parameters${MEASLES}.c ../../SIR.h
	$(CC) -c sampling_Parameters${MEASLES}.c $(CXXFLAGS) 

tempRate_SIR_PNAS_Version.o: ../tempRate_SIR_PNAS_Version.c ../../SIR.h
	$(CC) -c ../tempRate_SIR_PNAS_Version.c $(CXXFLAGS) 

init_SIR.o : ../../init_SIR.c ../../SIR.h
	$(CC) -c ../../init_SIR.c $(CXXFLAGS) 

IO_SIR.o : ../../IO_SIR.c ../../SIR.h
	$(CC) -c ../../IO_SIR.c $(CXXFLAGS) 

func_SIR.o:	../../func_SIR.c ../../SIR.h
	$(CC) -c ../../func_SIR.c $(CXXFLAGS)

func_Open_SIR.o: ../../PwSp_Param_Space/func_Open_SIR.c ../../SIR.h
	$(CC) -c ../../PwSp_Param_Space/func_Open_SIR.c $(CXXFLAGS) 

func_PwSp_Analytic.o: ../../PwSp_Param_Space/func_PwSp_Analytic.c ../../SIR.h
	$(CC) -c ../../PwSp_Param_Space/func_PwSp_Analytic.c $(CXXFLAGS) 

step_SIR.o:	../../step_SIR.c ../../SIR.h
	$(CC) -c ../../step_SIR.c $(CXXFLAGS) 

evolving_Parameters.o:	../../evolving_Parameters.c ../../SIR.h
	$(CC) -c ../../evolving_Parameters.c $(CXXFLAGS) 

ODEqs.o : ${ODE} ${HEAD.H} 
	${CC} -c ${ODE} $(CXXFLAGS)

ODE_definition.o : ../../TempEvo_Det/ODE_definition.c ../../SIR.h
	${CC} -c ../../TempEvo_Det/ODE_definition.c $(CXXFLAGS)

clean:  
	rm -f *.o SIR *~ core pFr*.dat inf*.dat sus*.dat          


