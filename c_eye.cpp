//
//
// File: c_eye.cpp
//
// Main file in solving the cylindrical eye problem. Initializes the 
// displacement functions R and W and solves for their values iteratively.
//
// @author Roland Sanford <ras9841@rit.edu>
// 
// // // // // // // // // // // // // // // // // // // // // // // // //

//
// Includes
//

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>

//
// Macros
//

#define NUM_ITER 5000
#define dr (R_EDGE/N)
#define dz (DEPTH/M)
#define OMEGA 1.41
#define OMEGA_LOW 1.38

//
// Globals
//

int M, N;
double E;		    // dynes/cm^2,	Young's modulus of eye	
double SIGMA;		//     	    	Poisson's ration of CL
double R_EDGE;	    // cm, 		    radius of undeformed CL
double DEPTH;       // cm,		    depth of eye
double DELTA;       // cm,          convergence condition
double *P;          // dynes

//
// Functions
//

void print_disp(double **function);
void write_csv(double **function, char ch);
void read_config(FILE *file, char *filename);
double function_average(double **W);
double r(int i, int j);

// Main functino in the cylindrical solution.
// Populates R and W with zeros as an initial
// guess and then runs the Gauss Seidel method.
//
// Preconditions:
//		none
// Postconditions:
//		R and W equilibrium values calculated 
//		for each point (i,j) on the grid.
int main(int argc, char *argv[]){
    if (argc != 2) {
        printf("Usage: ./clp config_file.txt\n");
        return EXIT_SUCCESS;
    }
    FILE *pFile;
    pFile = fopen(argv[1], "r");
    if(!pFile){
        printf("Error opening file.\n");
        return EXIT_FAILURE;
    }

    read_config(pFile, argv[1]); 
    
    // Populate R and W
	double **R = new double*[M+1];
	double **W = new double*[M+1];
    double **W_old = new double*[M+1];
    
	for (size_t i = 0; i < M + 1; i++){
		R[i] = new double[N+1];
		W[i] = new double[N+1];
		W_old[i] = new double[N+1];
        for (size_t j = 0; j < N+1; j++){
			R[i][j] = 0.00;
			W[i][j] = 0.00;
		}
	}
 
	// Setup timer
	std::clock_t start;
	double duration;
	start = std::clock();
	
    double curr_diff = 100.0, max_diff, old, gs, diff;
    size_t count = 0;
    
    while (curr_diff > DELTA){
		max_diff = 0;
        
        for (size_t i = 0; i < M+1; i++){
            memcpy(W_old[i], W[i], sizeof(double)*(N+1));         
        }

        // (0,0) (lower left corner
        R[0][0] = 0;
        
        old = W[0][0];
		gs = W[0][1];
        W[0][0] = old + OMEGA*(gs-old);
        
       
		// (0,M) (top left corner)
		R[M][0] = 0;
		
        old = W[M][0];
        gs = W[M][1]; 
        W[M][0] = old + OMEGA_LOW*(gs-old);
		
        // (N,0) (lower right corner) 
		R[0][N] = (
                    (1-SIGMA)*R[0][N-1]/dr - SIGMA*(W[1][N]-W[0][N])/dz
                  )
                  /( (1-SIGMA)/dr + SIGMA/r(0,N) ); 
		W[0][N] = W[0][N-1];

		// (N,M) (top right corner)
        old = R[M][N];
        gs = 4*R[M-1][N]/3 - R[M-2][N]/3 
                - dz*(3*W[M][N]-4*W[M][N-1]+W[M][N-2])/(3*dr); 
	    R[M][N] = old + OMEGA*(gs-old);
        if ( (diff = std::abs(old-R[M][N])) > max_diff ){
            max_diff = diff;
        }
        
        old = W[M][N];
        gs = 4*W[M-1][N]/3 - W[M-2][N]/3 
                - ( 2*dz*SIGMA/( 3*(1-SIGMA) ) )*(
                (3*R[M][N]-4*R[M][N-1]+R[M][N-2])/(2*dr)+R[M][N]/r(M,N))
                - 2*dz*P[N]*(1+SIGMA)*(1-2*SIGMA)/(3*(1-SIGMA)*E);
        W[M][N] = old + OMEGA*(gs-old);
        
 
		// i = 0 (lower bound w/o corners)
		for (int j = 1; j < N; j++){
			old = R[0][j];
            gs = R[1][j] + dz/(2*dr)*( W[0][j+1]-W[0][j-1] );
            R[0][j] = old + .8*(gs-old);
            if ( (diff = std::abs(old-R[0][j])) > max_diff ){
                max_diff = diff;
            }

            old = W[0][j];
            gs = W[1][j] + (SIGMA*dz / (1 - SIGMA)) *((R[0][j+1] - R[0][j-1]) 
                     / (2*dr) + R[0][j] / r(0, j));    
            W[0][j] = old + 1.1*(gs-old);
        }

		// i = M (top bound w/o corners)
		for (int j = 1; j < N; j++){
            old = R[M][j];
			gs = R[M-1][j] - dz *(W[M][j+1] - W[M][j-1]) / (2 * dr);
            R[M][j] = old + .8*(gs - old);
            if ( (diff = std::abs(old-R[M][j])) > max_diff ){
                 max_diff = diff;
            }

            old = W[M][j];
            gs = W[M-1][j] 
				- dz*(1 + SIGMA)*(1 - 2*SIGMA)*P[j]/((1 - SIGMA)*E) 
				- ( dz*SIGMA / (1 - SIGMA) )*
				(
					(R[M][j+1] - R[M][j-1])/(2*dr)
					+ R[M][j]/r(M,j) 
				);
            W[M][j] = old + .8*(gs - old);
		}
		
		// j = 0 (left bound w/o corners)
		for (int i = 1; i < M; i++){
			R[i][0] = 0;
			
            old = W[i][0];
            gs = W[i][1]; 
		    W[i][0] = old + OMEGA*(gs-old);  
        }

		// j = N (right boundary w/o corners)
		for (int i = 1; i < M; i++){
			old = R[i][N];
            gs = 
				( 
				(1-SIGMA)*(4*R[i][N-1]-R[i][N-2])/(2*dr)  
				- (SIGMA)/(2*dz)*( W[i+1][N] - W[i-1][N] )
				)
			 	/ ( 3*(1-SIGMA)/(2*dr) + SIGMA/r(i,N) );
			R[i][N] = old + OMEGA*(gs-old);
            if ( (diff = std::abs(old-R[i][N])) > max_diff ){
                 max_diff = diff;
            }
            
            old = W[i][N];
            gs = 4*W[i][N-1]/3 - W[i][N-2]/3 - (dr/(3*dz))*(R[i+1][N]- R[i-1][N]);
		    W[i][N] = old + OMEGA_LOW*(gs-old);
        }
		
        // Inside points
		for (int i = 1; i < M; i++){
			for (int j = 1; j < N; j++){
				old = R[i][j];
                gs  = (
					  2*(1-SIGMA)/(dr*dr)*(R[i][j+1]+R[i][j-1])
					+ 2*(1-SIGMA)/(2*dr*r(i,j))*(R[i][j+1]-R[i][j-1])  
					+ (1-2*SIGMA)*(R[i+1][j]+R[i-1][j])/( dz*dz )
					+ (W[i+1][j+1]-W[i+1][j-1]-W[i-1][j+1]+W[i-1][j-1])/( 4*dr*dz )
					) / 
					( 
					  4*(1-SIGMA)/( dr*dr )
					+ 2*(1-2*SIGMA)/( dz*dz )
					+ 1/( r(i,j)*r(i,j) )
					);
				R[i][j] = old + 1.86*(gs-old);
                if ( (diff = std::abs(old-R[i][j])) > max_diff ){
                    max_diff = diff;
                }

                old = W[i][j];
                gs  = 
					(
					  (W[i][j+1]+W[i][j-1])/(dr*dr)
					+ (W[i][j+1]-W[i][j-1])/(2*dr*r(i,j))
					+ 2*(1-SIGMA)/(dz*dz*(1-2*SIGMA))*(W[i+1][j]+W[i-1][j])
					+ ( 
					    (R[i+1][j]-R[i-1][j])/(2*dz*r(i,j)) 
					  + (R[i+1][j+1]-R[i+1][j-1]-R[i-1][j+1]+R[i-1][j-1])/(4*dr*dz) 
					  ) / (1-2*SIGMA) 
					)
					/ ( 2/(dr*dr)+4*(1-SIGMA)/(dz*dz*(1-2*SIGMA)) );
			    W[i][j] = old + 1.75*(gs-old);
            }
        }
        
        double av = function_average(W);
        for(size_t i = 0; i<M+1; i++){
            for(size_t j=0; j<N+1; j++){
                W[i][j] -= av;
                if ( (diff = std::abs(W_old[i][j]-W[i][j])) > max_diff ){ 
                    max_diff = diff;
                }
            }
        }

        //printf("Max diff: %e\n", max_diff);

        curr_diff = max_diff;
		count++;
	}	
	
	duration = (std::clock() - start)/(double)CLOCKS_PER_SEC;
    printf("Number of iterations: %zu\n", count);
    std::cout<<"Time:  "<< duration << "s. " <<"\n\n\n";

	write_csv(R, 'R');
	write_csv(W, 'W');
	
    #ifdef OCTAVE	
    std::system("octave display.m");
	#endif
	
	for (size_t i = 0; i < M + 1; i++){
		delete [] R[i];
		delete [] W[i];
	}	
	
	delete [] R;
	delete [] W;	
	delete [] P;
    return 0;
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
	for (size_t i = M; i < -1; i--){
		for (size_t j = 0; j < N+1; j++){
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
	for (size_t i = M; i < -1; i--){
		for (size_t j = 0; j < N+1; j++){
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

void read_config(FILE *pFile, char *filename){
    char buff[50], buff2[50];
    printf("\n\n\nUsing config file %s ...\n\n", filename);
    pFile = fopen(filename, "r");
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%d", &M);
    N = M;
    printf("%s %dx%d\n", buff, M, N);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &E);
    printf("%s %s %e\n", buff, buff2, E); 
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &SIGMA);
    printf("%s %g\n", buff, buff2, SIGMA);   
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &R_EDGE);
    printf("%s %s %g\n", buff, buff2, R_EDGE);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%s", &buff2);
    fscanf(pFile, "%lf", &DEPTH);
    printf("%s %s %g\n", buff, buff2, DEPTH);
    
    fscanf(pFile, "%s", &buff);
    fscanf(pFile, "%lf", &DELTA);
    printf("%s %e\n", buff, DELTA);  
    
    // Pressure
    P = new double[N+1];
    fscanf(pFile, "%s", &buff);
    for (size_t i = 0; i<N+1; i++){
        fscanf(pFile, "%lf", &P[i]);
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
