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



double Anl_tilde(int n ,int l){
    /*
    Function that computes the A_nl from Eq.16 in Lowing+16 paper

    Parameters:
    ----------
    n, l

    aeturn:
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


void sum_angular_prod(double * All_phi_mS, double * All_phi_mT, double * All_phi_mST, int n_points,\
                      double *r, double *theta, double *phi, double *M, \
                      int n, int l, int m, int n_prime, int l_prime,\
                      int m_prime){

    double phi_nlm, phi_nlm_prime;
    int i;
    *All_phi_mS = 0;
    *All_phi_mST = 0;
    *All_phi_mT = 0;

    for(i=0;i<=n_points;i++){
      phi_nlm = phi_nlm_f(r[i], theta[i], n, l, m);
      phi_nlm_prime = phi_nlm_f(r[i], theta[i], n_prime, l_prime, m_prime);
      *All_phi_mS += phi_nlm*phi_nlm_prime*cos(m*phi[i])*cos(m_prime*phi[i])*M[i]*M[i];
      *All_phi_mT += phi_nlm*phi_nlm_prime*sin(m*phi[i])*sin(m_prime*phi[i])*M[i]*M[i];
      *All_phi_mST += phi_nlm*phi_nlm_prime*sin(m*phi[i])*cos(m_prime*phi[i])*M[i]*M[i];
    }
}


void cov_matrix(int n_points, double *r , double *theta , double *phi,\
                double *M, int nmax, int lmax, char *out_filename){

    int n, l, m;
    double S_tilde[nmax+1][lmax+1][lmax+1];
    double T_tilde[nmax+1][lmax+1][lmax+1];
    double ST_tilde[nmax+1][lmax+1][lmax+1];

    #pragma omp parallel for private(l, m)
    for(n=0;n<=nmax;n++){
      printf("computing coefficient n= %d \n", n);
      for(l=0;l<=lmax;l++){
        double A_nl;
        double All_phi_nlm_mix_S;
        double All_phi_nlm_mix_T;
        double All_phi_nlm_mix_ST;
        A_nl = Anl_tilde(n,l);

        sum_angular_prod(&All_phi_nlm_mix_S, &All_phi_nlm_mix_T,\
                         &All_phi_nlm_mix_ST, n_points, r, theta, phi, M, n, l, 0,\
                         n, l, 0);
        S_tilde[n][l][0] = A_nl*A_nl*All_phi_nlm_mix_S;
        T_tilde[n][l][0] = A_nl*A_nl*All_phi_nlm_mix_T;
        ST_tilde[n][l][0] = A_nl*A_nl*All_phi_nlm_mix_ST;
        for(m=1;m<=l;m++){

           sum_angular_prod(&All_phi_nlm_mix_S, &All_phi_nlm_mix_T,\
                            &All_phi_nlm_mix_ST, n_points, r, theta, phi, M, n, l, m,\
                            n, l, m);

           S_tilde[n][l][m] = 4*A_nl*A_nl*All_phi_nlm_mix_S;
           T_tilde[n][l][m] = 4*A_nl*A_nl*All_phi_nlm_mix_T;
           ST_tilde[n][l][m] = 4*A_nl*A_nl*All_phi_nlm_mix_ST;
              
            
          
        }
      }
    }
   printf("Done computing covariance matrix and now writing the results \n");
   write_cov(out_filename,  nmax+1, lmax+1, S_tilde, T_tilde, ST_tilde);
   //write_mass(out_mass, M, n_points);

}


void coefficients(int n_points, double *r , double *theta , double *phi, double *M, int nmax, int lmax, char *out_filename){

    //
    //basic function that computes the S_nlm coefficient using the
    //HO basis function. Following Lowing et al. (2011).


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
   printf("Done computing coefficients and now writing the results \n");
   write_coeff(out_filename,  nmax+1, lmax+1, S, T);
}


