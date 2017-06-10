#ifndef TERMO_H
#define TERMO_H
double energia(int *lattice, int n, int j, int b, int seg);
double magnetizacion(int *lattice, int n);
double func_resp(double *vec_termo, int ns, int n, double *termo_e_m);
double correlacion(double *vec_termo, int nsamp, int n_vec_corr, double *vec_corr, int ns);
#endif
