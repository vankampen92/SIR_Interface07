#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*	                                                                      */
#/*                            Power Spectrum Analysis                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
C_DIR = /home/david/C_pro/
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip 
CXXFLAGS+=-I${C_DIR}C_NumRecip -g #-Wall  
OBJECTS = stat.o argumentControl.o report_PS.o main_PowerSPC.o IO_PS.o 
CC=gcc  

all: PWS

PWS:  	$(OBJECTS)   
	$(CC) -o PWS $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}/stat.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}stat.c $(CXXFLAGS)

main_PowerSPC.o: main_PowerSPC.c ../pwsp.h
	$(CC) -c main_PowerSPC.c $(CXXFLAGS)

IO_PS.o:../IO_PS.c ../pwsp.h
	${CC} -c ../IO_PS.c ${CXXFLAGS} 

argumentControl.o: argumentControl.c ../pwsp.h
	${CC} -c argumentControl.c $(CXXFLAGS)

report_PS.o : report_PS.c ../pwsp.h
	$(CC) -c report_PS.c $(CXXFLAGS) 

clean:  
	rm -f *.o PWS *~ core         


