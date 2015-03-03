//
// File: cylindrical_eye.cpp
//
// Main file in solving the cylindrical eye problem. Initializes the 
// displacement functions R and W and solves for their values by calling
// gs_method at each point in the computational grid.
//
// @author Roland Sanford <ras9841@rit.edu>
// 
// // // // // // // // // // // // // // // // // // // // // // // // //

//
// Includes
//

#include <iostream>
#include "gs_method.h"

//
// Globals
//

const int M = 100;
const int N = 100;
double E = 2.45;
double SIGMA = .4;


//
// Macros
//

#define NUM_ITER 1

//
// Functions
//

// Main functino in the cylindrical solution.
// Populates R and W with zeros as an initial
// guess and then runs the Gaus Seidel method.
//
// Preconditions:
//
// Postconditions:
//		R and W equilibrium values calculated 
//		for each point (i,j) on the grid.
int main(){

	// Populate R and W
	std::vector<std::vector<double>> R, W;

	R.resize(M, std::vector<double>(N, 0.0));
	W.resize(M, std::vector<double>(N, 0.0));


	int count = 0;
	while (count < NUM_ITER){
		for (int i = 0; i < M; i++){
			for (int j = 0; j < N; j++){
				gs_method(i, j, R, W);
			}
		}
		count++;
	}

	char a;
	std::cin >> a;
	return 0;
}