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
ODE_01 = ${C_DIR}/Differential_Equations/ODEqs_01.c 
PARAM = ${C_DIR}Stochastics/SIR/PwSp_Param_Space/ParamSpace_Numeric/
LDFLAGS+=-lm -lNumRec -lNRUTL -L${C_DIR}C_NumRecip -L${C_DIR}
CXXFLAGS+=-I${C_DIR}C_NumRecip -I${C_DIR} -g #-Wall 
OBJECTS = stat.o argumentControl.o report_SIR.o main_SIR.o func_SIR.o func_Open_SIR.o func_PwSp_Analytic.o init_SIR.o IO_eps.o eps_WGW_image.o z_vertical_scale.o evolving_Parameters.o peak_To_trough_Ratio.o ODE_definition.o ODEqs.o ODEqs_01.o det_seasonal_amplif.o sto_seasonal_amplif.o the_Largest_Maxima.o fig_TEX.o append_Fig_TEX.o model_SIR.o tempRate_SIR.o step_SIR.o
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

model_SIR.o: ${C_DIR}/Stochastics/SIR/TempEvo/model_SIR.c ../../../SIR.h
	$(CC) -c ${C_DIR}/Stochastics/SIR/TempEvo/model_SIR.c $(CXXFLAGS)

tempRate_SIR.o: ${C_DIR}/Stochastics/SIR/tempRate_SIR.c ../../../SIR.h
	$(CC) -c ${C_DIR}/Stochastics/SIR/tempRate_SIR.c $(CXXFLAGS) 

step_SIR.o:	${C_DIR}/Stochastics/SIR/step_SIR.c ../../../SIR.h
	$(CC) -c ${C_DIR}/Stochastics/SIR/step_SIR.c $(CXXFLAGS) 

det_seasonal_amplif.o: ${PARAM}det_seasonal_amplif.c ../../../SIR.h ../../../SIR_Analytic.h
	$(CC) -c ${PARAM}det_seasonal_amplif.c $(CXXFLAGS)

sto_seasonal_amplif.o: ${PARAM}sto_seasonal_amplif.c ../../../SIR.h ../../../SIR_Analytic.h
	$(CC) -c ${PARAM}sto_seasonal_amplif.c $(CXXFLAGS)

the_Largest_Maxima.o: ${PARAM}the_Largest_Maxima.c ../../../SIR.h ../../../SIR_Analytic.h
	$(CC) -c ${PARAM}the_Largest_Maxima.c $(CXXFLAGS)

argumentControl.o: argumentControl.c ../../../SIR.h
	${CC} -c argumentControl.c $(CXXFLAGS)

report_SIR.o: report_SIR.c ../../../SIR.h
	${CC} -c report_SIR.c $(CXXFLAGS)

func_SIR.o:	../../../func_SIR.c ../../../SIR.h
	$(CC) -c ../../../func_SIR.c $(CXXFLAGS) 

func_Open_SIR.o: ../../func_Open_SIR.c ../../../SIR.h ../../../SIR_Analytic_General.h
	$(CC) -c ../../func_Open_SIR.c $(CXXFLAGS) 

func_PwSp_Analytic.o: ../../func_PwSp_Analytic.c ../../../SIR.h ../../../SIR_Analytic_General.h
	$(CC) -c ../../func_PwSp_Analytic.c $(CXXFLAGS) 

init_SIR.o : ../../../init_SIR.c ../../../SIR.h
	$(CC) -c ../../../init_SIR.c $(CXXFLAGS) 

ODEqs.o : ${ODE} ${HEAD.H} 
	${CC} -c ${ODE} $(CXXFLAGS)

ODEqs_01.o : ${ODE_01} ${HEAD.H} 
	${CC} -c ${ODE_01} $(CXXFLAGS)

ODE_definition.o : ${C_DIR}Stochastics/SIR/TempEvo_Det/ODE_definition.c ${C_DIR}/Stochastics/SIR/SIR.h
	${CC} -c ${C_DIR}Stochastics/SIR/TempEvo_Det/ODE_definition.c $(CXXFLAGS)

evolving_Parameters.o:	${C_DIR}Stochastics/SIR/evolving_Parameters.c  ../../../SIR.h 
	$(CC) -c ${C_DIR}Stochastics/SIR/evolving_Parameters.c $(CXXFLAGS) 

peak_To_trough_Ratio.o : ${C_DIR}Stochastics/SIR/peak_To_trough_Ratio.c  ../../../SIR.h 
	$(CC) -c ${C_DIR}Stochastics/SIR/peak_To_trough_Ratio.c $(CXXFLAGS) 

clean:  
	rm -f *.o APWS *~ core         







