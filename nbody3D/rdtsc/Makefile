CC = mpicc 
CFLAGS = -O3 -Wall

all: measure_time


measure_time : measure_time.c
	${CC} ${CFLAGS} -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $^ -o $@  -lm -fopenmp


clean :	
	rm measure_time 
	
run : 
	mpirun ./measure_time
