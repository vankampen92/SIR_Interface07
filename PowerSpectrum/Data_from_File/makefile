#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*	                                                                      */
#/*                            Power Spectrum Analysis                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
HOME = /home/dalonso/
C_DIR = ${HOME}C_pro/
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -I${C_DIR} -g #-Wall  
OBJECTS = stat.o report_PW.o main_PW.o powerSpectrum.o argumentControl.o autoVar.o spectral_Density.o
CC=gcc  

all: PWS

PWS:  	$(OBJECTS)   
	$(CC) -o PWS $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}/stat.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}stat.c $(CXXFLAGS)

main_PW.o: main_PW.c ../pwsp.h
	$(CC) -c main_PW.c $(CXXFLAGS) 

argumentControl.o: argumentControl.c ../pwsp.h
	${CC} -c argumentControl.c $(CXXFLAGS)

spectral_Density.o: ${C_DIR}C_NumRecip/spectral_Density.c 
	${CC} -c ${C_DIR}C_NumRecip/spectral_Density.c $(CXXFLAGS)

report_PW.o : report_PW.c ../pwsp.h
	$(CC) -c report_PW.c $(CXXFLAGS) 

autoVar.o : ../autoVar.c ../pwsp.h
	$(CC) -c ../autoVar.c $(CXXFLAGS) 

powerSpectrum.o : ../powerSpectrum.c ../pwsp.h
	$(CC) -c ../powerSpectrum.c $(CXXFLAGS) 

clean:  
	rm -f *.o PWS *~ core         

rmevo: 
	rm pFr*.dat


