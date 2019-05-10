gcc -Wall -I/usr/local/include  -c main.c io.c  covariance_matrix.c   -lm -fopenmp -O2 
gcc -Wall -L/usr/local/lib main.o io.o covariance_matrix.o  -lgsl -lgslcblas -lm -fopenmp -O2 
#./a.out > cov_mat_sph_1x1.txt
