# Makefile for compiling the RR code on Linux systems

CC = gcc
#CF = -O3 -Wall -Wextra -Wno-unused-parameter -Wuninitialized -Winit-self -pedantic -std=gnu11
CF = -O3 -Wall -fopenmp  

INCLUDES = -I/usr/bin/include -I/home/rpaviot/Cuba-4.2/ 
LIBRARIES = -L/usr/bin/lib -lgsl -lgslcblas -lm -L/home/rpaviot/Cuba-4.2/ -lcuba -fopenmp
#CUBA_OBJ = /home/rpaviot/Cuba-4.2/*.o

#FINALIZE COMPILE FLAGS
CF += $(OPTIONS) #-g


## FINALIZE 
RR_INC = $(INCLUDES)
RR_LIB =  $(LIBRARIES) 



RR: compute_RR.c
	$(CC) compute_RR_v4.c -o compute_kernel $(CF) $(RR_INC) $(RR_LIB)


clean:
	rm -f compute_kernel */*~ *~

