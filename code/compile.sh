gcc -Wall -I/usr/local/include -c io.c covariance_matrix.c main.c   -lm -fopenmp -O2 -w
gcc -Wall -L/usr/local/lib io.o covariance_matrix.o main.o  -lgsl -lgslcblas -lm -fopenmp -O2 -w
