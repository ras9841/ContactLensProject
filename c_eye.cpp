//
// File: cylindrical_eye.cpp
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
#include <cmath>
#include <vector>

//
// Globals
//

const int M = 5;
const int N = 5;
const double E = std::pow(10,6);			// dynes/cm^s
const double SIGMA = .5;					//     Poisson's ration of CL
const double R_EDGE = .7;					// cm, radius of undeformed CL
const double H = 3 * std::pow(10,-4);		// cm, PLTF thickness
const double TAU = 4 * std::pow(10 ,-3);	// cm, thickness of undeformed CL

//
// Macros
//

#define NUM_ITER 5

#define dr 1.0/M
#define dz 1.0/N
#define r(i,j) j


//
// Functions
//

// Pressure function. Used to calculate the pressure due the a contact lens
// on the top boundary.
//
// Preconditions:
// 	none
// Postconditions:
// 	pressure at the point (i,j) calculated and returned.
double P(int i, int j){
	return 2;
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
	double W[M + 1][N + 1];

	// Initialize R and W
	for (size_t i = 0; i < M + 1; i++){
		for (size_t j = 0; j < N; j++){
			R[i][j] = 0.00;
			W[i][j] = 0.00;
		}
	}

	char a;
	for (size_t i = 0; i < M + 1; i++){
		for (size_t j = 0; j < N; j++){
			std::cout << R[i][j] << ", ";
		}
		std::cout << "\n";
	}
	std::cin >> a;
	
	int count = 0;
	while (count < NUM_ITER){
		// (0,0)
		R[0][0] = 0;
		W[0][0] = ((dr*dr) * (dz)) / ((dr*dr)*(1 - SIGMA) - 4 * dz*(1 - 2 * SIGMA))
			* (-4 / ((dr*dr))*(1 - 2 * SIGMA)*W[1][0] - (1 - SIGMA)
			/ (dz)*(3 * W[2][0] - 4 * W[1][0]));

		// (0,N)
		R[0][N] = 0;
		W[0][N] = W[0][N-1]-(2 * SIGMA*dz)/(1 - SIGMA) * R[1][N]/dr 
			- (dz*(1+SIGMA)*(1 - 2 * SIGMA)*P(0, N)) / ((1 - SIGMA)*E);

		// (M,0)
		R[M][0] = R[M][1] - dz / dr*W[M - 1][0];
		W[M][0] = 0;

		// (M,N)
		R[M][N] = -(((SIGMA - 1) / dr - (1 - SIGMA) / dz)*W[M][N-1] - ((SIGMA - 1)
			* (SIGMA - 1) / (SIGMA*dr) + SIGMA / dr)*R[M-1][N] + (1 + SIGMA)*(1 - 
			2 * SIGMA)*P(M, N) * 1 / E) / (((SIGMA - 1)*(SIGMA - 1) + (SIGMA*SIGMA)) 
			/ (SIGMA*dr) + 1 / r(M, N));
		W[M][N] = (-SIGMA*(1 / dr + 1 / (r(M, N)))*R[M][N] + (1 - SIGMA) / dz
			* W[M][N-1] + SIGMA / dr*R[M-1][N] - (1 + SIGMA)
			* (1 - 2 * SIGMA)*P(M, N) / E) / ((1 - SIGMA) / dz);

		
		// i = 0 (lower bound w/o corners)
		for (int j = 1; j < N; j++){
			R[0][j] = r(0, j)*(1 - SIGMA) / SIGMA*(R[0][j - 1] - R[0][j + 1]) 
				/ (2 * dr) - r(0, j)*(W[1][j] - W[0][j]) / dz;
			W[0][j] = W[1][j] + (SIGMA*dz / (1 - SIGMA)) *((R[0][j + 1] - R[0][j - 1]) 
				/ (2 * dr) + R[0][j] / r(0, j));
		}

		// i = M (top bound w/o corners)
		for (int j = 1; j < N; j++){
			R[M][j] = R[M-1][j] - dz*(W[M][j+1] - W[M][j-1]) / (2 * dr);
			W[M][j] = W[M-1][j] - dz*(1 + SIGMA)*(1 - 2 * SIGMA)*P(M, j) 
				/ ((1 - SIGMA)*E) - (dz*SIGMA / (1 - SIGMA))*((R[M][j+1] - R[M][j-1]) 
				/ (2 * dr)+R[M][j]/r(M,j));
		}

		// j = 0 (left bound w/o corners)
		for (int i = 1; i < M; i++){
			R[i][0] = 0;
			W[i][0] = (1/(dr*dr)*(2*W[i][1]-W[i][2])-((1-SIGMA)/(dz*dz*(1-2*SIGMA)))
				* (W[i+1][0]+W[i-1][0]))/(1/(dr*dr)-2*(1-SIGMA)/(dz*dz*(1-2*SIGMA)));
		}


		// j = N (right boundary w/o corners)
		for (int i = 1; i < M; i++){
			R[i][N] = (R[i][N-1]/dr - (1-SIGMA)/(2*SIGMA*dz)*(W[i+1][N] - W[i-1][N]))
				/ (1/dr + 1/r(i,N));
			W[i][N] = W[i][N-1]- (dr/(2*dz))*(R[i+1][N]- R[i-1][N]);
		}
		
		// Inside points
		for (int i = 1; i < M; i++){
			for (int j = 1; j < N; j++){
				R[i][j] = 0;
				W[i][j] = 0;
			}
		}

		count++;

		char a;
		for (size_t i = M; i < -1; i--){
			for (size_t j = 0; j < N+1; j++){
				std::cout << R[i][j] << ", ";
			}
			std::cout << "\n";
		}
		std::cin >> a;

	}
	std::cin >> a;

	return 0;
}



