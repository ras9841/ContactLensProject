#ifndef _HELPER_
#define _HELPER_

#include <stdio.h>
#include "spline.h"

extern double dr, dz, E, SIGMA, R_EDGE, EYE_EDGE;
extern int M, N;

void usage();
void print_disp(double**function);
void write_output(double *P, double *T, double *BIG_R);
void write_csv(double **function, char ch);
double r(int i, int j);
void read_config(char *filenames[], double **P, double *DEPTH, double *DELTA,
                 std::vector<double> *data_R, std::vector<double> *data_Z,
                 std::vector<double> *lens_R, std::vector<double> *lens_Z,
                 std::vector<double> *tau_R, std::vector<double> *tau_Z);
double function_average(double **W);

#endif
