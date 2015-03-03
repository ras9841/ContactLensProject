//
// File: gs_method.c
//
// Implementation of the Gauss-Seidel method. Solves for Rij and Wij
// based on location in the computational grid.
//
// @author Roland Sanford <ras9841@rit.edu>
//
// // // // // // // // // // // // // // // // // // // // // // // // //

//
// Includes
//

#include "gs_method.h"
#include <cmath>

//
// Macros
//

#define dr 1
#define dz 1
#define r(i,j) std::sqrt(i**2 + j**2)

//
// Functions
//

// Computes the value of Rij and Wij at the given point (i,j).
//
// Preconditions:
//		0 <= i <= M
//		0 <= j <= N
//		R_i and W are initialized to zero for the first iteration.
// Postconditions:
//		Rij and Wij values updated.
void gs_method(int i, int j, std::vector<std::vector<double>>& R, std::vector<std::vector<double>>& W){
	if (i == M){ // right boundary
		if (j == N){ // top right 
			R[M][N] = 0;
			W[M][N] = 0;
		}
		else if (j == 0){ //bottom right
			R[M][0] = 0;
			W[M][0] = 0;
		}
		else{ // general right boundary

			R[M][j] = -((2 * dr*SIGMA) / (1 - SIGMA)*((W[M][j + 1] - W[M][j - 1]) / (2 * dz)) - 4 * R[M - 1][j] + R[M - 2][j]) 
				/ (3 + (2 * dr*SIGMA) / (1 - SIGMA));

			W[M][j] = (4 * W[M - 1][j] - W[M - 2][j] - (dr / dz)*(R[M][j + 1] - R[M][j - 1])) / 3;

		}
	}

	else if (i == 0){
		if (j == 0){ // bottom left

		}
		else if (j == N){ // top left

		}
		else { // general left boundary.


		}
	}
}