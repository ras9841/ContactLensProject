#ifndef _HELPER_
#define _HELPER_

#include <stdio.h>

extern double dr, dz, E, SIGMA, R_EDGE;

void get_pressure(double *P, double *f, double *g, double *TAU, 
                  double *T_CL, double *R_CL, int M, int N);
void print_disp(double**function, int M, int N);
void write_csv(double **function, char ch, int M, int N);
double r(int i, int j);
void read_config(char *filenames[], int *M, int *N, double **P, double *E, 
                 double *SIGMA, double *R_EDGE, double *DEPTH, double *DELTA,
double **f, double **g, double **TAU);
double function_average(double **W, int M, int N);

#endif
