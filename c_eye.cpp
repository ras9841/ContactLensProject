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
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>

//
// Macros
//

#define NUM_ITER 100000
#define TOLERANCE std::pow(10,-15)
#define dr (R_EDGE/N)
#define dz (DEPTH/M)
#define PI 3.14159265


//
// Globals
//

const int M = 500;
const int N = 500;
const double E = std::pow(10,6);		// dynes/cm^2,	Young's modulus of eye	
const double SIGMA = .4;			//     		Poisson's ration of CL
const double R_EDGE = 1;			// cm, 		radius of undeformed CL
const double DEPTH =  1;			// cm,		depth of eye
const double H = 3 * std::pow(10,-4);		// cm, 		PLTF thickness
const double TAU = 4 * std::pow(10 ,-3);	// cm, 		thickness of undeformed CL


//
// Functions
//

void print_disp(double **function);
void write_csv(double **function, char ch);
double r(int i, int j);
double P(int i, int j);

// Main functino in the cylindrical solution.
// Populates R and W with zeros as an initial
// guess and then runs the Gauss Seidel method.
//
// Preconditions:
//		none
// Postconditions:
//		R and W equilibrium values calculated 
//		for each point (i,j) on the grid.
int main(){

	// Populate R and W
	double **R = new double*[M+1];
	double **W = new double*[M+1];
       
	for (size_t i = 0; i < M + 1; i++){
		R[i] = new double[N+1];
		W[i] = new double[N+1];
		for (size_t j = 0; j < N+1; j++){
			R[i][j] = 0.00;
			W[i][j] = 0.00;
		}
	}
	
	// Setup timer
	std::clock_t start;
	double duration;
	start = std::clock();
	
    double curr_diff = 100.0, max_diff, old, diff;
	size_t count = 0;
	while (curr_diff > std::pow(10,-11)){
		max_diff = 0;
        // (0,0) (lower left corner
        R[0][0] = 0;
        old = W[0][0];
		W[0][0] = W[0][1];
        if ( (diff = std::abs(old-W[0][0])) > max_diff ){ 
            max_diff = diff;
        }

		// (0,M) (top left corner)
		R[M][0] = 0;
		old = W[M][0];
        W[M][0] = W[M][1]; 
        if ( (diff = std::abs(old-W[M][0])) > max_diff ){
            max_diff = diff;
        }
		
        // (N,0) (lower right corner) 
		R[0][N] = 0; 
		W[0][N] = 0;

		// (N,M) (top right corner)
        old = R[M][N];
        R[M][N] = 4*R[M-1][N]/3 - R[M-2][N]/3 
                - dz*(3*W[M][N]-4*W[M][N-1]+W[M][N-2])/(3*dr); 
		if ( (diff = std::abs(old-R[M][N])) > max_diff ){
            max_diff = diff;
        }
        
        old = W[M][N];
        W[M][N] = 4*W[M-1][N]/3 - W[M-2][N]/3 
                - ( 2*dz*SIGMA/( 3*(1-SIGMA) ) )*(
                (3*R[M][N]-4*R[M][N-1]+R[M][N-2])/(2*dr)+R[M][N]/r(M,N))
                - 2*dz*P(M,N)*(1+SIGMA)*(1-2*SIGMA)/(3*(1-SIGMA)*E);
       	if ( (diff = std::abs(old-W[M][N])) > max_diff ){
            max_diff = diff;
        }
        
 
		// i = 0 (lower bound w/o corners)
		for (int j = 1; j < N; j++){
			R[0][j] = 0;
			W[0][j] = 0;		
        }

		// i = M (top bound w/o corners)
		for (int j = 1; j < N; j++){
            old = R[M][j];
			R[M][j] = R[M-1][j] - dz *(W[M][j+1] - W[M][j-1]) / (2 * dr);
			if ( (diff = std::abs(old-R[M][j])) > max_diff ){
                 max_diff = diff;
            }

            old = W[M][j];
            W[M][j] = W[M-1][j] 
				- dz*(1 + SIGMA)*(1 - 2*SIGMA)*P(M, j)/((1 - SIGMA)*E) 
				- ( dz*SIGMA / (1 - SIGMA) )*
				(
					(R[M][j+1] - R[M][j-1])/(2*dr)
					+ R[M][j]/r(M,j) 
				);
	        if ( (diff = std::abs(old-W[M][j])) > max_diff ){
                 max_diff = diff;
            }

		}
		
		// j = 0 (left bound w/o corners)
		for (int i = 1; i < M; i++){
			R[i][0] = 0;
			
            old = W[i][0];
            W[i][0] = W[i][1]; 
		    if ( (diff = std::abs(old-W[i][0])) > max_diff ){
                 max_diff = diff;
            }

        }


		// j = N (right boundary w/o corners)
		for (int i = 1; i < M; i++){
			old = R[i][N];
            R[i][N] = 
				( 
				(1-SIGMA)*(4*R[i][N-1]-R[i][N-2])/(2*dr)  
				- (SIGMA)/(2*dz)*( W[i+1][N] - W[i-1][N] )
				)
			 	/ ( 3*(1-SIGMA)/(2*dr) + SIGMA/r(i,N) );
			if ( (diff = std::abs(old-R[i][N])) > max_diff ){
                 max_diff = diff;
            }
            
            old = W[i][N];
            W[i][N] = 4*W[i][N-1]/3 - W[i][N-2]/3 - (dr/(3*dz))*(R[i+1][N]- R[i-1][N]);
		    if ( (diff = std::abs(old-W[i][N])) > max_diff ){
                 max_diff = diff;
            }
            
        }
		
		// Inside points
		for (int i = 1; i < M; i++){
			for (int j = 1; j < N; j++){
				old = R[i][j];
                R[i][j] = (
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
				if ( (diff = std::abs(old-R[i][j])) > max_diff ){
                    max_diff = diff;
                }

                old = W[i][j];
                W[i][j] = 
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
			    if ( (diff = std::abs(old-W[i][j])) > max_diff ){
                    max_diff = diff;
                }
            }
		}
        std::cout << "Max: " << max_diff << "\n";

        curr_diff = max_diff;
		count++;
	}	
	
	duration = (std::clock() - start)/(double)CLOCKS_PER_SEC;
    printf("Number of iterations: %zu\n", count);
    std::cout<<"Time:  "<< duration << "s. " <<'\n';

	#ifdef OCTAVE
	write_csv(R, 'R');
	write_csv(W, 'W');
	std::system("octave display.m");
	#endif
	
	for (size_t i = 0; i < M + 1; i++){
		delete [] R[i];
		delete [] W[i];
	}	
	
	delete [] R;
	delete [] W;	

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


// Pressure function. Used to calculate the pressure due the a contact lens
// on the top boundary.
//
// Preconditions:
// 	none
// Postconditions:
// 	pressure at the point (i,j) calculated and returned.
double P(int i, int j){
	return 1;
	//return std::pow(r(i,j),2) - std::pow(R_EDGE,2);
	//return sin(2*PI*r(i,j)/R_EDGE); 
	//return std::exp(r(i,j));
	//return 1 - 2*r(i,j)*r(i,j);
	//return ((E*std::pow(TAU,3)*56*H)/(12*(1-SIGMA*SIGMA)*std::pow(R_EDGE,4)));
}


