#include "helper.h"
#include <iostream>
#include <cmath>
#include <fstream>

// print_disp() displays the specified displacement function. All printed 
// values are rounded to 15 decimal places. The functions are printed in
// the following orientation:
//
// 	[M][0]  |---------------| [M][N]
// 		    |---------------|
// 		    |---------------|
// 		    |---------------|
//		    |---------------|
//	[0][0]	|---------------| [0][N]
//
// Preconditions:
// 	the display function corresponding to the input character has 
// 		been initialized.
// Postconditions:
// 	function printed.
void print_disp(double **function, int M, int N){
	for (int i = M; i < -1; i--){
		for (int j = 0; j < N+1; j++){
			if (std::abs(function[i][j]) == 0){// 1*std::pow(10,-15)){
				std::cout << "0.00000e+01,\t";
			}
			else{
				std::cout << function[i][j] << ",\t";
			}
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";
}

// write_csv() writes the specified function in the same manner that
// print_disp() prints.
//
// Preconditions:
// 	function is not null.
// Postconditions:
// 	all files are closed.
void write_csv(double **function, char ch, int M, int N){
	char name[6];
	name[0] = ch;
	name[1] = '.';
	name[2] = 'c';
	name[3] = 's';
	name[4] = 'v';
	name[5] = '\0';
	
	std::ofstream file;
	file.open(name);
	for (int i = M; i > -1; i--){
		for (int j = 0; j < N+1; j++){
			if (std::abs(function[i][j]) == 0){// 1*std::pow(10,-15)){
				file << "0.00000e+01,\t";
			}
			else{
				file << function[i][j] << ",\t";
			}
		}
		file << "\n";
	}
	file.close();
}

// Radial function. Used to calculate the distance in cm from the center
// of the eye (r=0).
//
// Preconditions:
// 	none
// Postconditions:
// 	radial distance from the point (i,j) to the center 
// 	calculated and returned
double r(int i, int j){
	return j * dr;
}

void read_config(FILE *pFile, char *filename, int *M, int *N, double **P,
                 double *E, double *SIGMA, double *R_EDGE, double *DEPTH,
                 double *DELTA){
    char buff[50], buff2[50];
    printf("\n\n\nUsing config file %s ...\n\n", filename);
    pFile = fopen(filename, "r");
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%d", M);
    *N = *M;
    printf("%s %dx%d\n", buff, *M, *N);
   
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", E);
    printf("%s %s %e\n", buff, buff2, *E); 
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", SIGMA);
    printf("%s %g\n", buff, buff2, *SIGMA);   
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", R_EDGE);
    printf("%s %s %g\n", buff, buff2, *R_EDGE);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", DEPTH);
    printf("%s %s %g\n", buff, buff2, *DEPTH);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%lf", DELTA);
    printf("%s %e\n\n", buff, *DELTA);  
    
    // Pressure
    *P = new double[*N+1];
    fscanf(pFile, "%s", &buff);
    for (int i = 0; i<*N+1; i++){
        fscanf(pFile, "%lf", &((*P)[i]));
    }
}

double function_average(double **W, int M, int N){
    double avg = 0.0;
    for(int i=0; i<M+1; i++){    
        for(int j=0; j<N+1; j++){
            avg += W[i][j];
        }
    }
    avg = avg/(M*N);
    return avg;
}
