CC = gcc
CF = -O3 -Wall -Wextra -Wno-unused-parameter -Wuninitialized -Winit-self -pedantic -ffast-math -std=gnu11

INCLUDES = -I/home/rpaviot/these/package/Cuba-4.2/
#LIBRARIES = -L/data/home/mbreton/gsl-2.6/lib/ -lgsl -lgslcblas -lm -L/dec/users/mabreton/Cuba-4.2/ -lcuba -lm

LIBRARIES = -lgsl -lgslcblas -lm -L/home/rpaviot/these/package/Cuba-4.2/ -lcuba

# COMPILATION FLAGS
CF += $(OPTIONS) #-g

# FINALIZE 
RR_INC = $(INCLUDES)
RR_LIB = $(LIBRARIES) 

RR: compute_RR_v4.c
	$(CC) compute_RR_final.c -o compute_kernel_final $(CF) $(RR_INC) $(RR_LIB)

clean:
	rm -f compute_kernel  */*~ *~
