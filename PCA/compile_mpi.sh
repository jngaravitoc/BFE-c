mpicc -Wall -I/usr/local/include -c covariance_matrix_mpi.c -lm
mpicc -Wall -L/usr/local/lib covariance_matrix_mpi.o -lgsl -lgslcblas -lm
mpirun -n 16 ./a.out
