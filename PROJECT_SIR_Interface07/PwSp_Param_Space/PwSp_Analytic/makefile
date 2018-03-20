#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*	                                                                      */
#/*                            Power Spectrum Analysis                        */
#/*                                                                           */
#/*                             David Alonso, 2005 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
# In order to build a different main_SIR program, please use make SUF={whatever_the_suffix is} !!
HOME    =/home/dalonso/
C_DIR   =${HOME}C_pro/
MAIN    =main_SIR${SUF}
OVERALL =../ParamSpace_Numeric/
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -I${C_DIR} -g #-Wall 
OBJECTS = argumentControl.o report_SIR.o ${MAIN}.o func_SIR.o func_Open_SIR.o func_PwSp_Analytic.o init_SIR.o overall_Amplification.o qromb_Accuracy.o default_values.c stat.o 
CC=gcc  

all: PSD

PSD:  	$(OBJECTS)   
	$(CC) -o PSD $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}stat.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}stat.c $(CXXFLAGS)

${MAIN}.o: ${MAIN}.c ../../SIR.h ../../SIR_Analytic.h
	$(CC) -c ${MAIN}.c $(CXXFLAGS)

argumentControl.o: argumentControl.c ../../SIR.h
	${CC} -c argumentControl.c $(CXXFLAGS)

report_SIR.o: report_SIR.c ../../SIR.h
	${CC} -c report_SIR.c $(CXXFLAGS)

default_values.o: default_values.c ../../SIR.h
	${CC} -c default_values.c ${CXXFLAGS}

overall_Amplification.o: ${OVERALL}overall_Amplification.c ../../SIR.h ../../SIR_Analytic_General.h
	${CC} -c ${OVERALL}overall_Amplification.c $(CXXFLAGS)

qromb_Accuracy.o: ${OVERALL}qromb_Accuracy.c ../../SIR.h ../../SIR_Analytic_General.h
	${CC} -c ${OVERALL}qromb_Accuracy.c $(CXXFLAGS)

func_SIR.o:	../../func_SIR.c ../../SIR.h
	$(CC) -c ../../func_SIR.c $(CXXFLAGS) 

func_Open_SIR.o: ../func_Open_SIR.c ../../SIR.h ../../SIR_Analytic_General.h
	$(CC) -c ../func_Open_SIR.c $(CXXFLAGS) 

func_PwSp_Analytic.o: ../func_PwSp_Analytic.c ../../SIR.h ../../SIR_Analytic_General.h
	$(CC) -c ../func_PwSp_Analytic.c $(CXXFLAGS) 

init_SIR.o : ../../init_SIR.c ../../SIR.h
	$(CC) -c ../../init_SIR.c $(CXXFLAGS) 

clean:  
	rm -f *.o PSD *~ pFr*.dat core*.* core.* core         








