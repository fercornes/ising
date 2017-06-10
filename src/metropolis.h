#ifndef METROPOLIS_H
#define METROPOLIS_H
int metropolis(int *lattice, int n, int j, int b, int seg, double *vec_exp, double *energ, double *m);
int pick_site(int *lattice, int n, int *indices);
int flip(int *lattice, int n, int j, int b, int seg, int *indices, double *vec_exp, double *energ, double *m);
#endif
