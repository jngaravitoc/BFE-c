/*
Author: J. Nicolas Garavito-Camargo
08/19/2016
University of Arizona.

Code to compute the covariance matrix of the
coefficients of the basis expansions functions.
Bases on Weinberg 1996 Algorithm.

Requirements:
-------------
gsl

Usage:
-------
./a.out nmax lmax #ofparticles

To-Do:
------
1. Parallelize code.


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>

double Anl_tilde(int n, int l);
double phi_nl_f(double r, int n, int l);
double phi_nlm_f(double r, double theta, int n, int l, int m);
double sum_angular(int n_points, double *r, double *theta, double *phi, double *M, int n, int l, int m);
double sum_angular_prod(int n_points, double *r, double *theta, double *phi, double *M, int n, int l, int m, int n_prime, int l_prime, int m_prime);
void read_data(char *filename, int n_points, double *r, double *theta, double *phi, double *m, double r_s);
void cov_matrix(int n_points, double *r , double *theta , double *phi, double *m, int max, int lmax);
//void coefficients(int n_points, double *r , double *theta , double *phi, double *m, int max, int lmax);

int main(int argc, char **argv){
     double *r=NULL;
     double *theta=NULL;
     double *phi=NULL;
     double *M=NULL;

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

     char filename[100]="data/spherical_halo.txt";
     double r_s = 40.85;

     // ------------------------

     printf("%s \n", filename);

     /* Allocating memory for pointers */
     r = malloc(n_points*sizeof(long double));
     theta = malloc(n_points*sizeof(long double));
     phi = malloc(n_points*sizeof(long double));
     M = malloc(n_points*sizeof(long double));

     read_data(filename, n_points, r, theta, phi, M, r_s);
     cov_matrix(n_points, r, theta, phi, M, nmax, lmax);;
     return 0;
}


double Anl_tilde(int n ,int l){
    double K_nl, factor, A_nl;
    double gamma_factor;

    K_nl = 0.5*n*(n+4.*l+3.) + (l+1.)*(2.*l+1.);
    factor = pow(2,8.*l+6.) / (4.*M_PI*K_nl);
    gamma_factor = pow(gsl_sf_gamma(2.0*l+1.5),2)/gsl_sf_gamma(n+4.*l+3.);
    A_nl =-factor* gsl_sf_fact(n)*(n+2.*l+1.5)*gamma_factor;
    return A_nl;
}


double phi_nl_f(double r, int n, int l){
    double factor, s, C_n;
    factor = pow(r, l) * pow((1.+r), (-2.*l-1.)) * pow(4*M_PI, 0.5);
    s = (r-1)/(r+1);
    C_n = gsl_sf_gegenpoly_n(n,2*l+1.5, s);
    return -factor*C_n;
}


/* Function that computes the potential see Eq. in Lowing*/
double phi_nlm_f(double r, double theta ,int n, int l, int m){
    double Y_lm, phi_nl;
    Y_lm = gsl_sf_legendre_sphPlm(l, m, cos(theta));
    phi_nl = phi_nl_f(r, n, l);
    return phi_nl*Y_lm;
}


/* function that sums the angular terms over all the particles */
double sum_angular(int n_points, double *r, double *theta, double *phi, double *M, int n, int l, int m){
    double all_angular=0;
    double phi_nlm;
    int i;
    for(i=0;i<=n_points;i++){
    phi_nlm = phi_nlm_f(r[i], theta[i], n, l, m);
    //printf("%f \n", phi_nlm);
    all_angular += phi_nlm*cos(m*phi[i])*M[i];
    }
    return all_angular;
}


double sum_angular_prod(int n_points, double *r, double *theta, double *phi, double *M, int n, int l, int m, int n_prime, int l_prime, int m_prime){
    double all_angular=0;
    double phi_nlm, phi_nlm_prime;
    int i;
    for(i=0;i<=n_points;i++){
      phi_nlm = phi_nlm_f(r[i], theta[i], n, l, m);
      phi_nlm_prime = phi_nlm_f(r[i], theta[i], n_prime, l_prime, m_prime);
      all_angular += phi_nlm*phi_nlm_prime*cos(m*phi[i])*cos(m_prime*phi[i])*M[i];
    }
    return all_angular;
}
void cov_matrix(int n_points, double *r , double *theta , double *phi, double *M, int nmax, int lmax){

    int n, l, m, dm0, n_prime, l_prime, m_prime, dm0_prime;
    double A_nl, A_nl_prime;
    double All_phi_nlm, All_phi_nlm_prime, All_phi_nlm_mix;
    double S, S_prime, SS_prime, S_cov_mat;

    for(n=0;n<=nmax;n++){
      for(l=0;l<=lmax;l++){
        for(m=0;m<=l;m++){
          for(n_prime=0;n_prime<=nmax;n_prime++){
            for(l_prime=0;l_prime<=lmax;l_prime++){
              for(m_prime=0;m_prime<=l_prime;m_prime++){

                if(m==0){
                  dm0=1.0;
                }
                else{
                  dm0=0.0;
                }

                if(m_prime==0){
                  dm0_prime=1.0;
                }
                else{
                  dm0_prime=0.0;
                }

                A_nl = Anl_tilde(n,l);
                A_nl_prime = Anl_tilde(n_prime,l_prime);

                All_phi_nlm = sum_angular(n_points, r, theta, phi, M, n, l, m);
                All_phi_nlm_prime = sum_angular(n_points, r, theta, phi, M, n_prime, l_prime, m_prime);
                All_phi_nlm_mix = sum_angular_prod(n_points, r, theta, phi, M, n, l, m, n_prime, l_prime, m_prime);
                S = (2-dm0)*A_nl*All_phi_nlm;
                S_prime = (2-dm0_prime)*A_nl_prime*All_phi_nlm_prime;
                SS_prime = (2-dm0)*(2-dm0_prime)*A_nl*A_nl_prime*All_phi_nlm_mix;
                S_cov_mat = SS_prime - S*S_prime;

                printf("%f \n", S_cov_mat);
              }
            }
          }
        }
      }
    }
}

/*
void coefficients(int n_points, double *r , double *theta , double *phi, double *M, int nmax, int lmax){

    //
    //basic function that computes the S_nlm coefficient using the
    //HO basis function. Following Lowing et al. (2011).

    //TO DO:
    //------
    //Implement S as an array/matrix.
    //Implement T computation.

    //

    int n, l, m, dm0;
    double A_nl;
    double All_angular;
    double S;
    for(n=0;n<=nmax;n++){
        for(l=0;l<=lmax;l++){
            for(m=0;m<=l;m++){

            if(m==0){
            dm0=1.0;
            }
            else{
            dm0=0.0;
            }

            A_nl = Anl_tilde(n,l);
            All_angular = sum_angular(n_points, r, theta, phi, M, n, l, m);
            S = (2-dm0)*A_nl*All_angular;

            //intf("%f \t \n", S);
            }
        }
    }
}
*/

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
