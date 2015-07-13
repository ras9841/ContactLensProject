#include "helper.h"
#include "spline.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

void printdp(double *func){
    for (int i=0; i<N; i++){
        if (i%10 == 0){ printf("\n"); }
	    std::cout << func[i] << ",\t";
    }
}

void tst(tk::spline f){
    printf("F[0.5]:\t%lf\n",f(.5));
}

void get_pressure(double *P, double *f, double *g, double *TAU,
                  double *T, double *R_EYE, double *R_DISP, double *W){
    // Spline for f
    std::vector<double> X(N+1), Y(N+1), Y_new(N+1);
    for (int j=0; j<N+1; j++){
        X[j] = R_EYE[j];
        Y[j] = f[j];
        Y_new[j] = f[j] + W[j];
    }
   
    tk::spline f_init, f_new;
    f_init.set_points(X,Y);
    f_new.set_points(X,Y_new);

    // Forward Euler for R
    // j=0 case
    R_EYE[1] = R_EYE[0] + dr*(1+(T[1] - T[0])/(dr*(1+SIGMA))); 
    for (int j = 1; j < N; j++){         
       R_EYE[j+1] = R_EYE[j] + dr*
        (
        (T[j]+(1+SIGMA)*r(M,j)-SIGMA*R_EYE[j])*std::sqrt(
            ( 1+std::pow((g[j+1]-g[j+1])/(2*dr),2) )
            / ( 1+std::pow((f_new(j)-f_new(j-1))/(R_EYE[j]-R_EYE[j-1]),2) )  
            )/r(M,j)             
        ); 
    }

    // SOR for T
}



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
void print_disp(double **function){
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
void write_csv(double **function, char ch){
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

void file_error(char *filename){
    printf("Error opening file %s.\n", filename);
    exit(1);   
}

void read_config(char *filenames[], int *M, int *N, double **P, double *E, 
                 double *SIGMA, double *R_EDGE, double *DEPTH, double *DELTA,
                 double **g, double **TAU, std::vector<double> *data_R,
                 std::vector<double> *data_Z){
    // check files (DONT USE filenames[0])
    FILE *pFile, *pressureFile, *clFile, *eyeFile, *tauFile;
    pFile = fopen(filenames[1], "r");
    pressureFile = fopen(filenames[2], "r");
    clFile = fopen(filenames[3], "r");
    eyeFile = fopen(filenames[4], "r");
    tauFile = fopen(filenames[5], "r");

    if (!pFile) { file_error(filenames[1]); }
    if (!pressureFile) { file_error(filenames[2]); }
    if (!clFile) { file_error(filenames[3]); }
    if (!eyeFile) { file_error(filenames[4]); }
    if (!tauFile) { file_error(filenames[5]); }
 
    
    char buff[50], buff2[50];
    printf("\n\n\nUsing config file %s ...\n\n", filenames[1]);
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
    fscanf(pFile, "%lf", &EYE_EDGE);
    printf("%s %s %g\n", buff, buff2, EYE_EDGE);

    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", DEPTH);
    printf("%s %s %g\n", buff, buff2, *DEPTH);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%lf", DELTA);
    printf("%s %e\n\n", buff, *DELTA);  
    
    // Pressure
    *P = new double[*N+1];
    for (int i = 0; i<*N+1; i++){
        fscanf(pressureFile, "%lf", &((*P)[i]));
    }

    // Contact Lens Shape
    *g = new double[*N+1];
    for (int i = 0; i<*N+1; i++){
        fscanf(clFile, "%lf", &((*g)[i]));
    }
    
    //  Eye Shape
    //*f = new double[*N+1];
    data_R->resize(*N+1);
    data_Z->resize(*N+1);

    double r, z;
    for (int i = 0; i<*N+1; i++){
        fscanf(eyeFile, "%lf%*c %lf%*c", &r, &z);
        data_R->push_back(r);
        data_Z->push_back(z);
    }
 
    // Lens Thickness
    *TAU = new double[*N+1];
    for (int i = 0; i<*N+1; i++){
        fscanf(tauFile, "%lf", &((*TAU)[i]));
    } 

}


double function_average(double **W){
    double avg = 0.0;
    for(int i=0; i<M+1; i++){    
        for(int j=0; j<N+1; j++){
            avg += W[i][j];
        }
    }
    avg = avg/(M*N);
    return avg;
}
