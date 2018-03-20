#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/*                                  SIR MODEL                                */
#/*	                                                                      */
#/*                          (CONSTANT COMMUNITY SIZE)                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
HOME = /home/dalonso
ODE_DEF= ODE_definition
HEAD.H = ../SIR.h
C_DIR = ${HOME}/C_pro/
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -g -DDETERMINISTIC -I${C_DIR} #-Wall 
CC = gcc
ODE = ${C_DIR}/Differential_Equations/ODEqs.c
ODE_01 = ${C_DIR}/Differential_Equations/ODEqs_01.c
OBJECTS = main_Det.o ${ODE_DEF}.o report_SIR.o argumentControl.o ODEqs_01.o ODEqs.o func_SIR.o stat.o evolving_Parameters.o

INT : 	$(OBJECTS)   
	$(CC) -o INT $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

main_Det.o : main_Det.c ${HEAD.H} ${C_DIR}/stat.h
	 $(CC) -c main_Det.c $(CXXFLAGS) 

stat.o : ${C_DIR}/stat.c ${C_DIR}/stat.h ${HEAD.H}
	 $(CC) -c ${C_DIR}/stat.c $(CXXFLAGS)  

argumentControl.o: argumentControl.c ${HEAD.H}
	${CC} -c argumentControl.c $(CXXFLAGS)

report_SIR.o : report_SIR.c ${HEAD.H}
	$(CC) -c report_SIR.c $(CXXFLAGS) 

func_SIR.o:	../func_SIR.c ../SIR.h
	$(CC) -c ../func_SIR.c $(CXXFLAGS) 

evolving_Parameters.o:	../evolving_Parameters.c ../SIR.h
	$(CC) -c ../evolving_Parameters.c $(CXXFLAGS) 

ODEqs.o : ${ODE} ${HEAD.H} 
	${CC} -c ${ODE} $(CXXFLAGS)

ODEqs_01.o : ${ODE_01} ${HEAD.H} 
	${CC} -c ${ODE_01} $(CXXFLAGS)

${ODE_DEF}.o : ../${ODE_DEF}.c ${HEAD.H}
	${CC} -c ../${ODE_DEF}.c $(CXXFLAGS)

clean:  
	rm -f *.o INT *~ core    
