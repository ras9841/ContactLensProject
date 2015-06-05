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


//
// Globals
//

const int M = 5;
const int N = 5;
const double E = std::pow(10,6);		// dynes/cm^2,	Young's modulus of eye	
const double SIGMA = .4;			//     		Poisson's ration of CL
const double R_EDGE = .7;			// cm, 		radius of undeformed CL
const double DEPTH = .7;
const double H = 3 * std::pow(10,-4);		// cm, 		PLTF thickness
const double TAU = 4 * std::pow(10 ,-3);	// cm, 		thickness of undeformed CL


//
// Macros
//

#define NUM_ITER 1
#define r(i,j) (j*dr)
#define dr (R_EDGE/N)
#define dz (DEPTH/M)


//
// Functions
//

void print_disp(double function[][N+1]);

// Pressure function. Used to calculate the pressure due the a contact lens
// on the top boundary.
//
// Preconditions:
// 	none
// Postconditions:
// 	pressure at the point (i,j) calculated and returned.
double P(int i, int j){
	return 0; //.49-r(i,j)*r(i,j);
	//return ((E*std::pow(TAU,3)*56*H)/(12*(1-SIGMA*SIGMA)*std::pow(R_EDGE,4)));
}

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
	double R[M+1][N+1];
	double W[M+1][N+1];
	
	std::cout << "Initial Guess: \n";

	// Initialize R and W
	for (size_t i = 0; i < M + 1; i++){
		for (size_t j = 0; j < N+1; j++){
			R[i][j] = 0.00;
			W[i][j] = 0.00;
		}
	}
	
	printf("R:\n");
	print_disp(R);	
	
	printf("W:\n");
	print_disp(W);	
	
	printf("-----> Final Solution after (%d) iterations\n", NUM_ITER);

	std::cout << "\n";
	
	int count = 0;
	while (count < NUM_ITER){
		// (0,0) (lower left corner)
		R[0][0] = 0;
		W[0][0] = W[1][0]+2*SIGMA*dz/((1-SIGMA)*dr)*R[0][1];

		// (0,M) (top left corner)
		R[M][0] = 0;
		W[M][0] = W[M-1][0] - (dz/(1-SIGMA))
			* ( 2*SIGMA/dr*R[M][1] +(1+SIGMA)*(1-2*SIGMA)*P(M,0)/E ); 

		// (N,0) (lower right corner) 
		R[0][N] = R[1][N] - (dz/dr)*W[0][N-1];
		W[0][N] = 0;

		// (N,M) (top right corner)
		R[M][N] = ( SIGMA*R[M][N-1]/dr - (1-SIGMA)*(W[M][N]- W[M-1][N])/dz  
			- (1+SIGMA)*(1-2*SIGMA)*P(M,N)/E )
			/ (SIGMA*(1/dr+1/r(M,N)));  
		W[M][N] = W[M][N-1] - (dr/dz)*(R[M][N]-R[M-1][N]);

		
		// i = 0 (lower bound w/o corners)
		for (int j = 1; j < N; j++){
			R[0][j] = r(0, j)*(1 - SIGMA) / SIGMA*(R[0][j-1] - R[0][j+1]) 
				/ (2*dr) - r(0, j)*(W[1][j] - W[0][j]) / dz;
			W[0][j] = W[1][j] + (SIGMA*dz / (1 - SIGMA)) *((R[0][j+1] - R[0][j-1]) 
				/ (2*dr) + R[0][j] / r(0, j));
		}

		// i = M (top bound w/o corners)
		for (int j = 1; j < N; j++){
			R[M][j] = R[M-1][j] - dz *(W[M][j+1] - W[M][j-1]) / (2 * dr);
			W[M][j] = W[M-1][j] - dz*(1 + SIGMA)*(1 - 2*SIGMA)*P(M, j) 
			 	/ ((1 - SIGMA)*E) - (dz*SIGMA / (1 - SIGMA))*((R[M][j+1] - R[M][j-1]) 
				/ (2 * dr)+R[M][j]/r(M,j));
		}

		// j = 0 (left bound w/o corners)
		for (int i = 1; i < M; i++){
			R[i][0] = 0;
			W[i][0] =  
				( 2/(dr*dr)*W[i][1]  - ( (1-SIGMA)/(dz*dz*(1-2*SIGMA)) )  
				* (W[i+1][0]+W[i-1][0]) ) 
				/ ( 2/(dr*dr)-2*(1-SIGMA)/(dz*dz*(1-2*SIGMA)) );
		}


		// j = N (right boundary w/o corners)
		for (int i = 1; i < M; i++){
			R[i][N] =  
				( (1-SIGMA)*R[i][N-1]/dr  - (SIGMA)/(2*dz)*(W[i+1][N] - W[i-1][N]))
			 	/ ((1-SIGMA)/dr + SIGMA/r(i,N));
			W[i][N] = W[i][N-1] - (dr/(2*dz))*(R[i+1][N]- R[i-1][N]);
		}
		
		// Inside points
		for (int i = 1; i < M; i++){
			for (int j = 1; j < N; j++){
				R[i][j] = (
					  2*(1-SIGMA)*(R[i][j+1]+R[i][j-1])/( dr*dr  ) 
					+ (1-SIGMA)*(R[i][j+1]-R[i][j-1])/( r(i,j)*dr  )
					+ (1-2*SIGMA)*(R[i+1][j]+R[i-1][j])/( dz*dz )
					+ (W[i+1][j+1]-W[i+1][j-1]-W[i-1][j+1]+W[i-1][j-1])/( 4*dr*dz )
					) / 
					( 
					  4*(1-SIGMA)/( dr*dr )
					+ 2*(1-2*SIGMA)/( dz*dz )
					+ 1/( r(i,j)*r(i,j) )
					);
				W[i][j] = 
					(
					(W[i][j+1]+W[i][j-1])*(1/(dr*dr)+ 1/(r(i,j)*2*dr))
					+ 2*(1-SIGMA)/(dz*dz*(1-2*SIGMA))*(W[i+1][j]+W[i-1][j])
					+ ( (R[i+1][j]-R[i-1][j])/(2*dz*r(i,j)) + 
					  (R[i+1][j+1]-R[i+1][j-1]-R[i-1][j+1]+R[i-1][j-1])/(4*dr*dz) )
						/(1-2*SIGMA) 
					)
					/ (2/(dr*dr)+4*(1-SIGMA)/(dz*dz*(1-2*SIGMA)));
			}
		}

		count++;
	}	
	printf("R:\n");
	print_disp(R);
	printf("W:\n");
	print_disp(W);
	return 0;
}

// print_disp() displays the specified displacement function. All printed 
// values are rounded to 15 decimal places. The functions are printed in
// the following orientation:
//
// 	[M][0]  |---------------| [M][N]
// 		|---------------|
// 		|---------------|
// 		|---------------|
//		|---------------|
//	[0][0]	|---------------| [0][N]
//
// Preconditions:
// 	the display function corresponding to the input character has 
// 		been initialized.
// Postconditions:
// 	function printed.
void print_disp(double function[][N+1]){
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


