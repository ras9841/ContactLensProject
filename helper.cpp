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

void usage() {
    printf("Usage: ./clp config_file.txt pressure_file.txt ");
    printf("lens_shape.txt eye_shape.txt tau.txt\n");
    exit(0);
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

void write_output(double *P, double *T, double *R)
{
    std::ofstream p, t, r;
    p.open("P.csv");
    t.open("T.csv");
    r.open("R.csv");
    for (int i =0; i<N+1; i++)
    {
        p << P[i] << ",\n";
        t << T[i] << ",\n";
        r << R[i] << ",\n";
    }
    p.close();
    t.close();
    r.close();
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

void read_config(char *filenames[], double **P, double *DEPTH, double *DELTA,
                 std::vector<double> *data_R, std::vector<double> *data_Z,
                 std::vector<double> *lens_R, std::vector<double> *lens_Z,
                 std::vector<double> *tau_R, std::vector<double> *tau_Z){
    
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
    printf("\n\n\nConfig File: \t\t\t%s\n", filenames[1]);
    printf("Eye Shape Input: \t\t%s\n", filenames[4]);
    printf("Lens Shape Input: \t\t%s\n", filenames[3]);
    printf("Thinkness Profile Input: \t%s\n\n", filenames[5]); 

    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%d", &M);
    N = M;
    printf("%s \t\t\t\t%dx%d\n", buff, M, N);
   
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &E);
    printf("%s %s \t\t%e\n", buff, buff2, E); 
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &SIGMA);
    printf("%s \t\t\t%g\n", buff, buff2, SIGMA);   
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &R_EDGE);
    printf("%s %s \t\t%g\n", buff, buff2, R_EDGE);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &EYE_EDGE);
    printf("%s %s \t\t%g\n", buff, buff2, EYE_EDGE);

    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", DEPTH);
    printf("%s %s \t\t\t%g\n", buff, buff2, *DEPTH);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%lf", DELTA);
    printf("%s \t\t\t%e\n\n", buff, *DELTA);  
    
    // Pressure
    *P = new double[N+1];
    for (int i = 0; i<N+1; i++){
        fscanf(pressureFile, "%lf", &((*P)[i]));
    }

    double r, z;
    //  Eye Shape
    while (fscanf(eyeFile, "%lf%*c %lf%*c", &r, &z) != EOF)
    {
        data_R->push_back(r);
        data_Z->push_back(z);
    }
    
    //  Lens Shape
    while (fscanf(clFile, "%lf%*c %lf%*c", &r, &z) != EOF)
    {   
        lens_R->push_back(r);
        lens_Z->push_back(z);
    }
    
    //  Lens Thickness (Tau)
    while (fscanf(tauFile, "%lf%*c %lf%*c", &r, &z) != EOF)
    {
        tau_R->push_back(r);
        tau_Z->push_back(z);
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
