#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*	                                                                      */
#/*                            Power Spectrum Analysis                        */
#/*                                                                           */
#/*                             David Alonso, 2000 (c)                        */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
HOME = /home/dalonso
C_DIR = ${HOME}/C_pro/
ODE = ${C_DIR}/Differential_Equations/ODEqs.c 
ODE_01 = ${C_DIR}/Differential_Equations/ODEqs_01.c 
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -I${C_DIR} -g #-Wall 
OBJECTS = stat.o argumentControl.o report_SIR.o main_SIR.o func_SIR.o func_Open_SIR.o func_PwSp_Analytic.o init_SIR.o IO_eps.o eps_WGW_image.o z_vertical_scale.o fig_TEX.o append_Fig_TEX.o resilience.o 
CC=gcc  

all: APWS

APWS:  	$(OBJECTS)   
	$(CC) -o APWS $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)  

stat.o:	${C_DIR}stat.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}stat.c $(CXXFLAGS)

IO_eps.o: ${C_DIR}IO_eps.c ${C_DIR}stat.h
	  ${CC} -c ${C_DIR}IO_eps.c $(CXXFLAGS)

eps_WGW_image.o: ${C_DIR}eps_WGW_image.c ${C_DIR}stat.h
	  ${CC} -c ${C_DIR}eps_WGW_image.c $(CXXFLAGS)

z_vertical_scale.o: ${C_DIR}z_vertical_scale.c ${C_DIR}stat.h
	  ${CC} -c ${C_DIR}z_vertical_scale.c $(CXXFLAGS)

fig_TEX.o: ${C_DIR}fig_TEX.c ${C_DIR}stat.h
	${CC} -c ${C_DIR}fig_TEX.c $(CXXFLAGS)

append_Fig_TEX.o: append_Fig_TEX.c ../../../SIR.h ../../../SIR_Analytic.h
	$(CC) -c append_Fig_TEX.c $(CXXFLAGS)

main_SIR.o: main_SIR.c ../../../SIR.h ../../../SIR_Analytic.h
	$(CC) -c main_SIR.c $(CXXFLAGS)

argumentControl.o: argumentControl.c ../../../SIR.h
	${CC} -c argumentControl.c $(CXXFLAGS)

report_SIR.o: report_SIR.c ../../../SIR.h
	${CC} -c report_SIR.c $(CXXFLAGS)

func_SIR.o:	../../../func_SIR.c ../../../SIR.h
	$(CC) -c ../../../func_SIR.c $(CXXFLAGS) 

resilience.o:	../../resilience.c ../../../SIR.h
	$(CC) -c ../../resilience.c $(CXXFLAGS) 

func_Open_SIR.o: ../../func_Open_SIR.c ../../../SIR.h ../../../SIR_Analytic_General.h
	$(CC) -c ../../func_Open_SIR.c $(CXXFLAGS) 

func_PwSp_Analytic.o: ../../func_PwSp_Analytic.c ../../../SIR.h ../../../SIR_Analytic_General.h
	$(CC) -c ../../func_PwSp_Analytic.c $(CXXFLAGS) 

init_SIR.o : ../../../init_SIR.c ../../../SIR.h
	$(CC) -c ../../../init_SIR.c $(CXXFLAGS) 

clean:  
	rm -f *.o APWS *~ core         







