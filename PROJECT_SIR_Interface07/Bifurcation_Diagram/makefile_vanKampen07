#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
#/*                                                                              *\
#/*                                SIR / SEIR MODEL                              * |
#/*	                                                                         * |
#/*                              (BIFURCATION DIAGRAM)                           * |
#/*                                                                              * |
#/*                              David Alonso, 2000 (c)                          * | 
#/*                                                                              * |
#/*    Warning: ETA=_Eta, to calculate the bifurcation diagram with              * |  
#/*    the transformation:                                                       * |
#/*                                                                              * |
#/*                Beta = (1 - p) Beta_0                                         * |
#/*                Eta  = (1 - p) Eta_0                                          * |
#/*                                                                              */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
HOME=/home/dalonso
ODE_DEF=ODE_definition
ETA=_Eta
C_DIR=${HOME}/C_pro/
HEAD.H=../SIR.h
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -g -DDETERMINISTIC -I${C_DIR}#-Wall 
CC = gcc
ODE = ${C_DIR}Differential_Equations
ODE_01 = ${ODE}/ODEqs_Double.c
OBJECTS = main_Det.o ${ODE_DEF}.o report_SIR.o argumentControl.o ODEqs_Double.o func_SIR.o stat.o evolving_Parameters.o poincareMap${ETA}.o default_values.o ODE_definitionSEIR.o ODE_definitionSIR.o euler.o

POINC : $(OBJECTS)   
	$(CC) -o POINC $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

main_Det.o : main_Det.c ${HEAD.H} ${C_DIR}/stat.h
	 $(CC) -c main_Det.c $(CXXFLAGS) 

stat.o : ${C_DIR}/stat.c ${C_DIR}/stat.h ${HEAD.H}
	 $(CC) -c ${C_DIR}/stat.c $(CXXFLAGS)

argumentControl.o: argumentControl.c ${HEAD.H}
	${CC} -c argumentControl.c $(CXXFLAGS)

default_values.o: default_values.c ${HEAD.H}
	${CC} -c default_values.c $(CXXFLAGS)

report_SIR.o : report_SIR.c ${HEAD.H}
	$(CC) -c report_SIR.c $(CXXFLAGS)

func_SIR.o:	../func_SIR.c ../SIR.h
	$(CC) -c ../func_SIR.c $(CXXFLAGS)

evolving_Parameters.o:	../evolving_Parameters.c ../SIR.h
	$(CC) -c ../evolving_Parameters.c $(CXXFLAGS)

poincareMap${ETA}.o: poincareMap${ETA}.c ${HEAD.H}
	${CC} -c poincareMap${ETA}.c $(CXXFLAGS)

ODEqs_Double.o : ${ODE_01} ${HEAD.H}
	${CC} -c ${ODE_01} $(CXXFLAGS)

euler.o : ${ODE}/euler.c ${HEAD.H}
	${CC} -c ${ODE}/euler.c $(CXXFLAGS)

${ODE_DEF}.o : ../ODE_definitions/${ODE_DEF}.c ${HEAD.H}
	${CC} -c ../ODE_definitions/${ODE_DEF}.c $(CXXFLAGS)

ODE_definitionSEIR.o : ../ODE_definitions/ODE_definitionSEIR.c ${HEAD.H}
	${CC} -c ../ODE_definitions/ODE_definitionSEIR.c $(CXXFLAGS)

ODE_definitionSIR.o : ../ODE_definitions/ODE_definitionSIR.c ${HEAD.H}
	${CC} -c ../ODE_definitions/ODE_definitionSIR.c $(CXXFLAGS)

clean:  
	rm -f *.o POINC *~ core*    











