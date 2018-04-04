# ==========================================
# CUSTUMIZED GSL-based UTILITY LIBRARY make
# ==========================================
#
# ==========================================
# Gnu C make to build GSL-based UTILITY Library 
# The generated library is libDa_GSL_UTL.a 
# ==========================================

GSL_Stat_UTL=${HOME}/PROJECT_STOCHASTIC_VTD_MODELS/Library/GSL_RANDOM/
EXTRA_INCLUDE_PATH= -I${GSL_Stat_UTL}
CFLAGS=-g -DHAVE_INLINE ${EXTRA_INCLUDE_PATH}
###CFLAGS= -Wall -O21 
COMPILER=gcc

OBJECTS = GSL_Discrete_Sampling.o GSL_Init_Random_Seed.o GSL_stat.o fprintf_to_File_Matrix_gsl.o	

.c.o: 
	${COMPILER} $(CFLAGS) -c $<

GSL_UTL: ${OBJECTS} GSL_stat.h
	ar -q libDa_GSL_UTL.a ${OBJECTS}

clean:  
	rm -f *.o *~
