/*
Author: J. Nicolas Garavito-Camargo
08/19/2016
University of Arizona.

Code to compute the covariance matrix of the
coefficients of the basis expansions functions.
Bases on Weinberg 1996 Algorithm.

Requirements:
-------------
gnu/gsl
openmp

Usage:
-------
./a.out nmax lmax #ofparticles snapshot outputfile

To-Do:
------
1. Put comments and organize code!
2. Make # of particles not an argument.
3. Read gadget snapshots.

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
#include "covariance_matrix.h"

int main(int argc, char **argv){
     double *r=NULL;
     double *theta=NULL;
     double *phi=NULL;
     double *M=NULL;

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

     /* Allocating memory for pointers */
     r = malloc(n_points*sizeof(long double));
     theta = malloc(n_points*sizeof(long double));
     phi = malloc(n_points*sizeof(long double));
     M = malloc(n_points*sizeof(long double));
     
     start = clock();
     read_data(argv[4], n_points, r, theta, phi, M, r_s);
     end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
     printf("time to load the data:  %f \n", cpu_time_used);
     //coefficients(n_points, r, theta, phi, M, nmax, lmax, argv[5]);
     cov_matrix(n_points, r, theta, phi, M, nmax, lmax, argv[5]);
     return 0;
}


double Anl_tilde(int n ,int l){
    /*
    Function that computes the A_nl from Eq.16 in Lowing+16 paper

    Parameters:
    ----------
    n, l

    Return:
    ------- 

    A_nl
    */
    
    double K_nl, factor, A_nl;
    double gamma_factor;

    // Definition of K_nl
    K_nl = 0.5*n*(n+4.*l+3.) + (l+1.)*(2.*l+1.);

    factor = pow(2,8.*l+6.) / (4.*M_PI*K_nl);
    gamma_factor = pow(gsl_sf_gamma(2.0*l+1.5),2) / gsl_sf_gamma(n+4.*l+3.);
    A_nl =-factor * gsl_sf_fact(n) * (n+2.*l+1.5) * gamma_factor;

    return A_nl;
}


double phi_nl_f(double r, int n, int l){
    /*
    Function that computes Phi_nl using Eq.11 in Lowing+16 paper.

    Parameters:
    -----------
    r : Distances  
    n :
    l

    Returns:
    -------
    Phi_nl
    */
    double factor, s, C_n;
    factor = pow(r, l) * pow((1.+r), (-2.*l-1.)) * pow(4*M_PI, 0.5);
    s = (r-1)/(r+1);
    C_n = gsl_sf_gegenpoly_n(n,2*l+1.5,s);
    return -factor*C_n;
}


/* Function that computes the potential phi_nlm in Lowing+11*/
double phi_nlm_f(double r, double theta ,int n, int l, int m){
    double Y_lm, phi_nl;
    Y_lm = gsl_sf_legendre_sphPlm(l, m, cos(theta));
    phi_nl = phi_nl_f(r, n, l);
    return phi_nl*Y_lm;
}


void sum_angular(double * All_phi_S, double * All_phi_T, int n_points, double *r, double *theta, double *phi,\
                 double *M, int n, int l, int m){

    double phi_nlm;
    int i;
    *All_phi_S = 0;
    *All_phi_T = 0;

    for(i=0;i<=n_points;i++){
    phi_nlm = phi_nlm_f(r[i], theta[i], n, l, m);
    *All_phi_S += phi_nlm*cos(m*phi[i])*M[i];
    *All_phi_T += phi_nlm*sin(m*phi[i])*M[i];

    }
}


void sum_angular_prod(double * All_phi_mS, double * All_phi_mT, int n_points,\
                      double *r, double *theta, double *phi, double *M, \
                      int n, int l, int m, int n_prime, int l_prime,\
                      int m_prime){

    double phi_nlm, phi_nlm_prime;
    int i;
    *All_phi_mS = 0;
    *All_phi_mT = 0;

    for(i=0;i<=n_points;i++){
      phi_nlm = phi_nlm_f(r[i], theta[i], n, l, m);
      phi_nlm_prime = phi_nlm_f(r[i], theta[i], n_prime, l_prime, m_prime);
      *All_phi_mS += phi_nlm*phi_nlm_prime*cos(m*phi[i])*cos(m_prime*phi[i])*M[i]*M[i];
      *All_phi_mT += phi_nlm*phi_nlm_prime*sin(m*phi[i])*sin(m_prime*phi[i])*M[i]*M[i];
    }
}


void cov_matrix(int n_points, double *r , double *theta , double *phi,\
                double *M, int nmax, int lmax, char *out_filename){

    int n, l, m;
    double S_tilde[nmax+1][lmax+1][lmax+1];
    double T_tilde[nmax+1][lmax+1][lmax+1];

    #pragma omp parallel for private(l, m)
    for(n=0;n<=nmax;n++){
      for(l=0;l<=lmax;l++){
        double A_nl;
        double All_phi_nlm_mix_S;
        double All_phi_nlm_mix_T;
        A_nl = Anl_tilde(n,l);

        sum_angular_prod(&All_phi_nlm_mix_S, &All_phi_nlm_mix_T,\
                         n_points, r, theta, phi, M, n, l, 0,\
                         n, l, 0);
        S_tilde[n][l][0] = A_nl*A_nl*All_phi_nlm_mix_S;
        T_tilde[n][l][0] = A_nl*A_nl*All_phi_nlm_mix_T;
        for(m=1;m<=l;m++){

           sum_angular_prod(&All_phi_nlm_mix_S, &All_phi_nlm_mix_T,\
                            n_points, r, theta, phi, M, n, l, m,\
                            n, l, m);

           S_tilde[n][l][m] = 4*A_nl*A_nl*All_phi_nlm_mix_S;
           T_tilde[n][l][m] = 4*A_nl*A_nl*All_phi_nlm_mix_T;
              
            
          
        }
      }
    }

   write_data(out_filename,  nmax+1, lmax+1, S_tilde, T_tilde);

   //for(n=0;n<=nmax;n++){
   //  for(l=0;l<=lmax;l++){
   //    for(m=0;m<=l;m++){
   //       printf("%8.4e \t %8.4e \t %d \t %d \t %d \n", S_tilde[n][l][m], T_tilde[n][l][m], n, l, m); 
   //     }
   //   }
   // }

}


void coefficients(int n_points, double *r , double *theta , double *phi, double *M, int nmax, int lmax, char *out_filename){

    //
    //basic function that computes the S_nlm coefficient using the
    //HO basis function. Following Lowing et al. (2011).

    //TO DO:
    //------
    //Implement S as an array/matrix.
    //Implement T computation.

    //

    int n, l, m;
    double S[nmax+1][lmax+1][lmax+1];
    double T[nmax+1][lmax+1][lmax+1];

    

    #pragma omp parallel for private(l, m)
    for(n=0;n<=nmax;n++){
       for(l=0;l<=lmax;l++){

            double A_nl;
            double All_phi_nlm_S, All_phi_nlm_T;

            A_nl = Anl_tilde(n,l);
            sum_angular(&All_phi_nlm_S, &All_phi_nlm_T, n_points, r, theta, phi, M, n, l, 0); 
            S[n][l][0] = A_nl*All_phi_nlm_S;  
            T[n][l][0] = A_nl*All_phi_nlm_T;  
            //printf("%f \t %f \n", S[n][l][0], T[n][l][0]); 

            for(m=1;m<=l; m++){
                       
            sum_angular(&All_phi_nlm_S, &All_phi_nlm_T, n_points, r, theta, phi, M, n, l, m);
            S[n][l][m] = 2.0*A_nl*All_phi_nlm_S;
            T[n][l][m] = 2.0*A_nl*All_phi_nlm_T;
        
            }
        }
    }
    //ite_data(out_filename, nmax, lmax, S , T); 
   write_data(out_filename,  nmax+1, lmax+1, S, T);
   //for(n=0;n<=nmax;n++){
   //  for(l=0;l<=lmax;l++){
   //    for(m=0;m<=l;m++){
   //       printf("%f \t %f \n", S[n][l][m], T[n][l][m]); 
  // }
   //}
   //}
}


void read_data(char *filename, int n_points, double *r, double *theta, \
               double *phi, double *M, double r_s){

    FILE *in;
    double X, Y, Z, m;
    int i;

    in = fopen(filename, "r");

    /* Checking if the file is opening*/
    if(!in){
        printf("Problem opening file %s \n", filename);
        exit(1);
    }

    for(i=0;i<=n_points;i++){
        fscanf(in, "%lf %lf %lf %lf \n", &X, &Y, &Z, &m);
        r[i] = pow(pow(X, 2) + pow(Y, 2) + pow(Z,2),0.5)/r_s;
        theta[i] = acos(Z/(r[i]*r_s));
        phi[i] = atan2(Y, X);
        M[i] = m;
    }
    fclose(in);
}





void write_data(char *filename, int n_max, int l_max, double S[n_max][l_max][l_max], double T[n_max][l_max][l_max]){

    FILE *out;
    int n, l, m;
    out = fopen(filename, "w");

    /* Checking if the file is opening*/
    if(!out){
        printf("Problem opening file %s \n", filename);
        exit(1);
    }

    fprintf(out, "#S \t T \t  n \t l \t m \n");
    for(n=0;n<n_max;n++){
       for(l=0;l<l_max;l++){
           for(m=0;m<=l;m++){
        fprintf(out, "%8.4e \t  %8.4e \t %d \t %d \t %d \n", S[n][l][m], T[n][l][m], n, l, m);
        }
      }
    }
    fclose(out);
}
