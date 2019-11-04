#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.h"


void read_data(char *filename, int n_points, double *r, double *theta, \
               double *phi, double *M, double r_s){

    FILE *in;
    double X, Y, Z, m;
    double VX, VY, VZ;
    int i;

    in = fopen(filename, "r");

    /* Checking if the file is opening*/
    if(!in){
        printf("Problem opening file %s \n", filename);
        exit(1);
    }

    for(i=0;i<=n_points;i++){
        fscanf(in, "%lf %lf %lf %lf %lf %lf %lf \n", &X, &Y, &Z, &VX, &VY, &VZ, &m);
        r[i] = pow(pow(X, 2) + pow(Y, 2) + pow(Z,2),0.5)/r_s;
        theta[i] = acos(Z/(r[i]*r_s));
        phi[i] = atan2(Y, X);
        M[i] = m;
    }
    fclose(in);
}





void write_coeff(char *filename, int n_max, int l_max, double S[n_max][l_max][l_max], double T[n_max][l_max][l_max]){

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


void write_cov(char *filename, int n_max, int l_max, double S[n_max][l_max][l_max], double T[n_max][l_max][l_max], double ST[n_max][l_max][l_max]){

    FILE *out;
    int n, l, m;
    out = fopen(filename, "w");

    /* Checking if the file is opening*/
    if(!out){
        printf("Problem opening file %s \n", filename);
        exit(1);
    }

    fprintf(out, "#S \t T \t  ST \t n \t l \t m \n");
    for(n=0;n<n_max;n++){
       for(l=0;l<l_max;l++){
           for(m=0;m<=l;m++){
        fprintf(out, "%8.4e \t %8.4e \t  %8.4e \t %d \t %d \t %d \n", S[n][l][m], T[n][l][m], ST[n][l][m], n, l, m);
        }
      }
    }
    fclose(out);
}

void rand_sampling(int n_points, double *r_rand, double *theta_rand, double *phi_rand, double *M_rand, double *r, double *theta, double *phi, double *M){
    srand(time(NULL));   
    int n;
    int rn;
    printf("performing random sampling \n");
    for(n=0;n<=n_points;n++){
        rn = rand()%n_points;      
        r_rand[n] = r[rn];
        theta_rand[n] = theta[rn];
        phi_rand[n] = phi[rn];
        M_rand[n] = M[rn];
    }
}

void read_data_sample(char *filename, int n_points, double *r, double *theta, \
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

    int n;
    int rn;

    for(n=0;n<=n_points;n++){
        rn = rand()%n_points;      
        //printf("%d \n", rn);
        //printf("%f \n", r[rn]);
        printf("%f \t %f \t %f \n", r[rn], theta[rn], phi[rn]);

    }
}
