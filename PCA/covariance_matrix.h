double Anl_tilde(int n, int l);                                                                                                                                          
double phi_nl_f(double r, int n, int l);                                                                                                                                 
double phi_nlm_f(double r, double theta, int n, int l, int m);                                                                                                           
void sum_angular(double * All_phi_S, double * All_phi_T, int n_points, double *r, double *theta, double *phi, double *M, int n, int l, int m);                                                                 
void sum_angular_prod(double * All_phi_S, double * All_phi_T, int n_points,\                                                                                             
                      double *r, double *theta, double *phi, double *M, \
                     int n, int l, int m, int n_prime, int l_prime, int m_prime); 
                                                                                                                                                                         
void read_data(char *filename, int n_points, double *r, double *theta,\                                                                                                  
               double *phi, double *m, double r_s);

void write_data(char *filename, int n_points, double *S);

void cov_matrix(int n_points, double *r , double *theta , double *phi,\                                                                                                  
                    double *m, int max, int lmax);   
void coefficients(int n_points, double *r , double *theta, double *phi, double *m, int max, int lmax);                                                                   

