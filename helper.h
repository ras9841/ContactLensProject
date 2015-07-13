#ifndef _HELPER_
#define _HELPER_

#include <stdio.h>
#include "spline.h"

extern double dr, dz, E, SIGMA, R_EDGE, EYE_EDGE;
extern int M, N;

void usage();
void print_disp(double**function);
void write_csv(double **function, char ch);
double r(int i, int j);
void read_config(char *filenames[], int *M, int *N, double **P, double *E, 
                 double *SIGMA, double *R_EDGE, double *DEPTH, double *DELTA,
                 double **g, double **TAU, std::vector<double> *data_R, 
                 std::vector<double> *data_Z);
double function_average(double **W);

#endif
