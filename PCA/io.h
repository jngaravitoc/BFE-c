void read_data(char *filename, int n_points, double *r, double *theta, double *phi, double *m, double r_s);

void write_data(char *filename, int n_max, int l_max, double S[n_max][l_max][l_max], double T[n_max][l_max][l_max]);

void rand_sampling(int n_points, double *r_rand, double *theta_rand, double *phi_rand, double *M_rand, double *r, double *theta, double *phi, double *M);
