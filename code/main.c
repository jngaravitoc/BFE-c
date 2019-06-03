/* Code written by Nico Garavito-Camargo
 * University of Arizona, 2017-2019
 * github: jngaravitoc/SCF_tools
 */


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>

#include "io.h"
#include "covariance_matrix.h"



int main(int argc, char **argv){
     double *r=NULL;;
     double *theta=NULL;
     double *phi=NULL;
     double *M=NULL;

     double *r_rand=NULL;
     double *theta_rand=NULL;
     double *phi_rand=NULL;
     double *M_rand=NULL;

     clock_t start, end;
     double cpu_time_used;


     if(!argv[1]){
     printf("Usage: ./a.out nmax lmax #ofparticles \n");
     }


     /* Global variables */
     int nmax= atof(argv[1]);
     int lmax= atof(argv[2]);
     // TO-DO: The code should recognize the number of particles in
     // the file.
     int n_points=atoi(argv[3]);
     //char filename;
     //filename = atol(argv[4]);
     printf("reading data %s /n", argv[4]);
     printf("WTF /n");
     //char filename[100];
     //filename = argv[4]; //"../data/spherical_halo.txt";
     double r_s = 40.85;
     // ------------------------


     int n, n_sampling;
     n_sampling = pow(n_points,0.5);
       

     /* Allocating memory for pointers */
     r = malloc(n_points*sizeof(long double));
     theta = malloc(n_points*sizeof(long double));
     phi = malloc(n_points*sizeof(long double));
     M = malloc(n_points*sizeof(long double));


     r_rand = malloc(n_points*sizeof(long double));
     theta_rand = malloc(n_points*sizeof(long double));
     phi_rand = malloc(n_points*sizeof(long double));
     M_rand = malloc(n_points*sizeof(long double));


     start = clock();
     read_data(argv[4], n_points, r, theta, phi, M, r_s);
     end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
     printf("time to load the data:  %f \n", cpu_time_used);
     //rand_sampling(n_points, r_rand, theta_rand, phi_rand, M_rand, r, theta, phi, M);
     coefficients(n_points, r, theta, phi, M, nmax, lmax, argv[5]);
     //coefficients(n_points, r_rand, theta_rand, phi_rand, M_rand, nmax, lmax, "coeff_rand.txt");

     cov_matrix(n_points, r, theta, phi, M, nmax, lmax, "cov_matrix_test.txt");

     //char buffer[93];
    
     /*
     for(n=0;n<=300;n++){
     printf("Computing covariance matrix using random sampling this is trial %d out of %d \n", n, n_sampling);
     snprintf(buffer, sizeof(char) * 93, "/home/u9/jngaravitoc/codes/BFE/covariance_bootstrap/cov_matrix_rand_mwlmc5_b1_1E6_%5i.txt", n);
     rand_sampling(n_points, r_rand, theta_rand, phi_rand, M_rand, r, theta, phi, M);
     cov_matrix(n_points, r_rand, theta_rand, phi_rand, M_rand, nmax, lmax, buffer);
     }*/
     return 0;
}
