CC = gcc
OPTIMIZE = -O2 -Wall -w 
EXEC = BFE
MATH_LIBS = -lm 
GSL_LIBS = -lgsl -lgslcblas
OPM_LIBS = -fopenmp
LDFLAGS = -L/usr/local/include
CFLAGS = -I/usr/local/lib
LFLAGS = $(GSL_LIBS) $(MATH_LIBS) $(OPM_LIBS)

INCL = allvars.h io.h covariance_matrix.h
OBJS = main.o io.o covariance_matrix.o



#############
#build rules#
#############

%.o: %.c $(INCL)
	$(CC) -c $(CFLAGS) -o $@ $< $(LFLAGS) 


$(EXEC): $(OBJS)
	$(CC) -c $(OPTIMIZE) $(LDFLAGS) $(LFLAGS) $(OBJS) $(LIBS) -o $(EXEC)  

.PHONY: clean

clean: 
	rm -f $(OBJS) 

