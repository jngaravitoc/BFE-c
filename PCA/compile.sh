gcc -Wall -I/usr/local/include -c covariance_matrix.c -lm
gcc -Wall -L/usr/local/lib covariance_matrix.o -lgsl -lgslcblas -lm
#./a.out > cov_mat_sph_1x1.txt
