/* Code written by Nico Garavito-Camargo
 * University of Arizona, 2017-2019
 * github: jngaravitoc/BFE-c
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>

#include "io.h"
#include "covariance_matrix.h"

int main(int argc, char **argv){


     printf("This is BFE-c running \n");

     double *r=NULL;
     double *theta=NULL;
     double *phi=NULL;
     double *M=NULL;

     double *r_rand=NULL;
     double *theta_rand=NULL;
     double *phi_rand=NULL;
     double *M_rand=NULL;

     clock_t start, end;
     double cpu_time_used;


     /* Global variables */
     int nmax=atoi(argv[1]);
     int lmax=atoi(argv[2]);
     // TODO: The code should recognize the number of particles in
     // the file.
     int covariance=atoi(argv[3]);
     int n_points=atoi(argv[4]);
     int n_samples=atoi(argv[5]);

     char *in_path = argv[6];
     char *in_snap = argv[7];
     char *out_path = argv[8];
     char *out_snap = argv[9];


     int n_sampling = atoi(argv[10]);
     int init_snap = atoi(argv[11]);
     int final_snap = atoi(argv[12]);
     double r_s = atof(argv[13]);

     int n_snaps;

     
     //sampling = pow(n_points,0.5);
       
     
     /* Allocating memory for pointers */
     r = malloc(n_points*sizeof(long double));
     theta = malloc(n_points*sizeof(long double));
     phi = malloc(n_points*sizeof(long double));
     M = malloc(n_points*sizeof(long double));


     r_rand = malloc(n_points*sizeof(long double));
     theta_rand = malloc(n_points*sizeof(long double));
     phi_rand = malloc(n_points*sizeof(long double));
     M_rand = malloc(n_points*sizeof(long double));

     for(n_snaps=init_snap;n_snaps<=final_snap;n_snaps++){
     
     start = clock();

     char buffer_str_in[10];
     char in_file[500];
     char mass_file[500];
     char ext_file[50];
     strcpy(in_file, in_path);
     strcat(in_file, in_snap);
     snprintf(ext_file, sizeof(char)*50, "_%03i.txt", n_snaps);
     strcat(in_file, ext_file);
     printf("reading data %s \n", in_file);

     read_data(in_file, n_points, r, theta, phi, M, r_s);
     end = clock();

     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
     printf("time to load the snapshot:  %f \n", cpu_time_used);

     char buffer_str_coeff[60], buffer_str_cov[60];
     char out_coeff[500], out_covmat[500];
    

     if(n_samples>0){     
          
     printf("Sampling particles \n");
     printf("Computing covariance matrix using random sampling this is trial %d  \n", n_samples);

     snprintf(buffer_str_coeff, sizeof(char) * 60, "_coeff_sample_%04i_snap_%04i.txt", n_samples, n_snaps);
     snprintf(buffer_str_cov, sizeof(char) * 60, "_covmat_sample_%04i_snap_%04i.txt", n_samples, n_snaps);
     
     strcpy(out_coeff, out_path);
     strcpy(out_covmat, out_path);
     strcat(out_coeff, out_snap);
     strcat(out_covmat, out_snap);
     strcat(out_coeff, buffer_str_coeff);
     strcat(out_covmat, buffer_str_cov);

     rand_sampling(n_sampling, r_rand, theta_rand, phi_rand, M_rand, r, theta, phi, M);
     printf("Done sampling halo \n");
     printf("Computing coefficients \n");
     coefficients(n_sampling, r_rand, theta_rand, phi_rand, M_rand, nmax, lmax, out_coeff, out_covmat, covariance);
     
     }     
     
     else{
     int n;
     n=0;

     snprintf(buffer_str_coeff, sizeof(char) * 60, "_coeff_sample_%04i_snap_%04i.txt", n, n_snaps);
     snprintf(buffer_str_cov, sizeof(char) * 60, "_covmat_sample_%04i_snap_%04i.txt", n, n_snaps);
     
     strcpy(out_coeff, out_path);
     strcpy(out_covmat, out_path);
     strcat(out_coeff, out_snap);
     strcat(out_covmat, out_snap);
     strcat(out_coeff, buffer_str_coeff);
     strcat(out_covmat, buffer_str_cov);
     
  
     printf("Computing coefficients \n");
     coefficients(n_points, r, theta, phi, M, nmax, lmax, out_coeff, out_covmat, covariance);
     
     printf("Computing covariance matrix \n");
     //cov_matrix(n_points, r, theta, phi, M, nmax, lmax, out_covmat);
     printf("Done computing covariance matrix \n");

     }
     } 
     
     return 0;
}
