#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*	                                                                      */
#/*                            Power Spectrum Analysis                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
HOME = /homes/dalonso/
C_DIR = ${HOME}C_pro/
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -I${C_DIR} -g #-Wall  
OBJECTS = stat.o argumentControl.o report_SIR.o main_SIR.o IO_SIR.o powerSpectrum.o model_SIR.o func_SIR.o func_Open_SIR.o init_SIR.o timeSeries.o timeSeries_VanKampen.o step_SIR.o tempRate_SIR_PNAS_Version.o autoCorr.o autoVar.o evolving_Parameters.o
CC=gcc  

all: APWS

APWS:  	$(OBJECTS)   
	$(CC) -o APWS $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}stat.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}stat.c $(CXXFLAGS)

main_SIR.o: main_SIR.c ../SIR.h
	$(CC) -c main_SIR.c $(CXXFLAGS)

argumentControl.o: argumentControl.c ../SIR.h
	${CC} -c argumentControl.c $(CXXFLAGS)

powerSpectrum.o: ${C_DIR}/PowerSpectrum/powerSpectrum.c ${C_DIR}/PowerSpectrum/pwsp.h
	${CC} -c ${C_DIR}/PowerSpectrum/powerSpectrum.c $(CXXFLAGS)

autoCorr.o: ${C_DIR}/PowerSpectrum/autoCorr.c ${C_DIR}/PowerSpectrum/pwsp.h
	${CC} -c ${C_DIR}/PowerSpectrum/autoCorr.c $(CXXFLAGS)

autoVar.o: ${C_DIR}/PowerSpectrum/autoVar.c ${C_DIR}/PowerSpectrum/pwsp.h
	${CC} -c ${C_DIR}/PowerSpectrum/autoVar.c $(CXXFLAGS)

model_SIR.o:	${C_DIR}/Stochastics/SIR/PNAS_DushoffEtal/model_SIR.c ../SIR.h
	$(CC) -c ${C_DIR}/Stochastics/SIR/PNAS_DushoffEtal/model_SIR.c $(CXXFLAGS)

report_SIR.o : report_SIR.c ../SIR.h
	$(CC) -c report_SIR.c $(CXXFLAGS) 

init_SIR.o : ../init_SIR.c ../SIR.h
	$(CC) -c ../init_SIR.c $(CXXFLAGS) 

evolving_Parameters.o : ../evolving_Parameters.c  ../SIR.h
	$(CC) -c ../evolving_Parameters.c $(CXXFLAGS)

timeSeries_VanKampen.o : ../timeSeries_VanKampen.c ../SIR.h
	$(CC) -c ../timeSeries_VanKampen.c $(CXXFLAGS) 

timeSeries.o : ../timeSeries.c ../SIR.h
	$(CC) -c ../timeSeries.c $(CXXFLAGS) 

IO_SIR.o : ../IO_SIR.c ../SIR.h
	$(CC) -c ../IO_SIR.c $(CXXFLAGS) 

tempRate_SIR_PNAS_Version.o: ${C_DIR}/Stochastics/SIR/PNAS_DushoffEtal/tempRate_SIR_PNAS_Version.c ../SIR.h
	$(CC) -c ${C_DIR}/Stochastics/SIR/PNAS_DushoffEtal/tempRate_SIR_PNAS_Version.c $(CXXFLAGS) 

func_SIR.o:	../func_SIR.c ../SIR.h
	$(CC) -c ../func_SIR.c $(CXXFLAGS) 

func_Open_SIR.o: ../PwSp_SIR_General/func_Open_SIR.c ../SIR.h
	$(CC) -c ../PwSp_SIR_General/func_Open_SIR.c $(CXXFLAGS) 

step_SIR.o:	../step_SIR.c ../SIR.h
	$(CC) -c ../step_SIR.c $(CXXFLAGS) 

clean:  
	rm -f *.o APWS *~ pFr*.dat core         



