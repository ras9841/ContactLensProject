#ifndef _HELPER_
#define _HELPER_

#include <stdio.h>

extern double dr; 
extern double dz;

void print_disp(double**function, int M, int N);
void write_csv(double **function, char ch, int M, int N);
double r(int i, int j);
void read_config(FILE *pFile, char *filename, int *M, int *N, double **P,
                 double *E, double *SIGMA, double *R_EDGE, double *DEPTH,
                 double *DELTA);
double function_average(double **W, int M, int N);

#endif
