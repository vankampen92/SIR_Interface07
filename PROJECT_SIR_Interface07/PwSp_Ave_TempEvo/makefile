#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/*                                  SIR MODEL                                */
#/*                                                                           */
#/*                          (CONSTANT COMMUNITY SIZE)                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
C_DIR = /home/dalonso/C_pro
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}/C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}/C_NumRecip -I${C_DIR} -g #-Wall  
OBJECTS = IO_SIR.o stat.o argumentControl.o report_SIR.o model_SIR.o func_SIR.o main_SIR.o init_SIR.o step_SIR.o tempRate_SIR.o evolving_Parameters.o powerSpectrum.o func_Open_SIR.o autoVar.o spectral_Density.o
CC=gcc  

all: SIR_TE

SIR_TE: $(OBJECTS)   
	$(CC) -o SIR_TE $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}/stat.c ${C_DIR}/stat.h
	${CC} -c ${C_DIR}/stat.c $(CXXFLAGS)

main_SIR.o: main_SIR.c ../SIR.h
	$(CC) -c main_SIR.c $(CXXFLAGS) 

model_SIR.o: ../TempEvo/model_SIR.c ../SIR.h
	$(CC) -c ../TempEvo/model_SIR.c $(CXXFLAGS)

report_SIR.o : report_SIR.c ../SIR.h
	$(CC) -c report_SIR.c $(CXXFLAGS) 

init_SIR.o : ../init_SIR.c ../SIR.h
	$(CC) -c ../init_SIR.c $(CXXFLAGS) 

IO_SIR.o : ../IO_SIR.c ../SIR.h
	$(CC) -c ../IO_SIR.c $(CXXFLAGS) 

tempRate_SIR.o: ../tempRate_SIR.c ../SIR.h
	$(CC) -c ../tempRate_SIR.c $(CXXFLAGS) 

powerSpectrum.o: ${C_DIR}/PowerSpectrum/powerSpectrum.c ${C_DIR}/PowerSpectrum/pwsp.h
	${CC} -c ${C_DIR}/PowerSpectrum/powerSpectrum.c $(CXXFLAGS)

spectral_Density.o: ${C_DIR}/C_NumRecip/spectral_Density.c 
	${CC} -c ${C_DIR}/C_NumRecip/spectral_Density.c $(CXXFLAGS)

autoVar.o: ${C_DIR}/PowerSpectrum/autoVar.c ${C_DIR}/PowerSpectrum/pwsp.h
	${CC} -c ${C_DIR}/PowerSpectrum/autoVar.c $(CXXFLAGS)

func_Open_SIR.o: ../PwSp_Param_Space/func_Open_SIR.c ../SIR.h
	$(CC) -c ../PwSp_Param_Space/func_Open_SIR.c $(CXXFLAGS) 

func_SIR.o:	../func_SIR.c ../SIR.h
	$(CC) -c ../func_SIR.c $(CXXFLAGS) 

step_SIR.o:	../step_SIR.c ../SIR.h
	$(CC) -c ../step_SIR.c $(CXXFLAGS) 

evolving_Parameters.o:	../evolving_Parameters.c ../SIR.h
	$(CC) -c ../evolving_Parameters.c $(CXXFLAGS) 

argumentControl.o: argumentControl.c ../SIR.h
	${CC} -c argumentControl.c $(CXXFLAGS)

clean:  
	rm -f *.o SIR_TE *~ core          




