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
#include <cmath>
#include <vector>

//
// Globals
//

int M = 100;
int N = 100;
double E = 10^(6);				// dynes/cm^s
double SIGMA = .5;				//     Poisson's ration of CL
double R_EDGE = .7;				// cm, radius of undeformed CL
double H = 3 * 10 ^ (-4);		// cm, PLTF thickness
double TAU = 4 * 10 ^ (-3);		// cm, thickness of undeformed CL

//
// Macros
//

#define NUM_ITER 1

#define dr 1/M
#define dz 1/N
#define r(i,j) (i)*dr
#define P(i,j) (i) //((E*TAU^3*56*H)/(12*(1-SIGMA^2)*(R_EDGE)^4))

//
// Functions
//

// Main functino in the cylindrical solution.
// Populates R and W with zeros as an initial
// guess and then runs the Gaus Seidel method.
//
// Preconditions:
//		none
// Postconditions:
//		R and W equilibrium values calculated
//		for each point (i,j) on the grid.
int main(){

	// Populate R and W
	std::vector<std::vector<double>> R, W;

	R.resize(M+1, std::vector<double>(N+1, 0.0));
	W.resize(M+1, std::vector<double>(N+1, 0.0));


	int count = 0;
	while (count < NUM_ITER){
		// (0,0)
		R[0][0] = 0;
		W[0][0] = ((dr)^2*(dz))/((dr*dr)*(1-SIGMA)-4*dz*(1-2*SIGMA))*(-4/((dr)^2)*(1-2*SIGMA)*W[1][0]-(1-SIGMA)/(dz)*(3*W[2][0]-4*W[1][0]));

		// (0,N)
		R[0][N] = 0;
		W[0][N] = W[0][N - 1] -(2*SIGMA*dz)/(1-SIGMA)*R[1][N]/dr- (dz*(1+SIGMA)*(1-2*SIGMA)*P(0,N))/((1-SIGMA)*E);

		// (M,0)
		R[M][0] = R[M][1]-dz/dr*W[M-1][0];
		W[M][0] = 0;

		// (M,N)
		R[M][N] = -(((SIGMA - 1) / dr - (1 - SIGMA) / dz)*W[M][N - 1] - ((SIGMA - 1)*(SIGMA - 1) / (SIGMA*dr) + SIGMA /dr)*R[M-1][N]
			+ (1+SIGMA)*(1-2*SIGMA)*P(M,N)*1/E) / (((SIGMA - 1)*(SIGMA - 1) + (SIGMA*SIGMA)) / (SIGMA*dr) + 1 / r(M, N));
		W[M][N] = ( -SIGMA*(1/dr + 1/(r(M,N)))*R[M][N]+(1-SIGMA)/dz*W[M][N-1]+SIGMA/dr*R[M-1][N]-(1+SIGMA)*(1-2*SIGMA)*P(M,N)/E )/((1-SIGMA)/dz);


		// j = 0 (lower bound w/o corners)
		for (int i = 1; i < M; i++){
			R[i][0] = (1 / 3)*(dz/dr*(W[i+1][0]-W[i-1][0]) - R[i][2] + 4*R[i][1]);
			W[i][0] = (-1 / 3)*((2 * dz / (1 - SIGMA))*((1 - 2 * SIGMA) / (2 * dz))*(4 * W[i][1] - W[i][2])
				+ SIGMA*((R[i + 1][0] - R[i - 1][0]) / (2 * dr)+R[i][0]/r(i,0)+(4*W[i+1][0]-3*W[i][2])/(2*dz)));
		}

		// i = M (right boundary w/o corners)
		for (int j = 1; j < N; j++){
			R[M][j] = -((2 * dr*SIGMA) / (1 - SIGMA)*((W[M][j + 1] - W[M][j - 1]) / (2 * dz)) - 4 * R[M - 1][j] + R[M - 2][j])
				/ (3 + (2 * dr*SIGMA) / (1 - SIGMA));

			W[M][j] = (4 * W[M - 1][j] - W[M - 2][j] - (dr / dz)*(R[M][j + 1] - R[M][j - 1])) / 3;
		}

		// Inside points
		for (int i = 1; i < M; i++){
			for (int j = 1; i < N; j++){
				double T = (R[i + 1][j] + R[i - 1][j]) / (dr ^ 2) * (1 + 1 / (1 - 2 * SIGMA)) - (R[i][j + 1] + R[i][j - 1]) / (dz ^ 2);
				double S = (1 / r(i, j))*( (R[i][j+1]-R[i-1][j])/(2*dr) + (W[i+1][j+1]-W[i-1][j+1]+W[i-1][j-1])/(4*dr*dz) );

				R[i][j] = ( (R[i+1][j]-R[i-1][j])/(2*dr*r(i,j)) + S/(1-2*SIGMA) +
					 T + (R[i][j+1]+R[i][j-1])/(dz^2) )
					/ ( (2 / dr)*(1+1/(1-2*SIGMA))+2/(dz^2)+1/r(i,j) );

				double H = -1/r(i,j)*(W[i+1][j] - W[i-1][j])/dr - (R[i][j+1]+R[i][j-1])/((1-2*SIGMA)*r(i,j)*2*dr)+
					(R[i+1][j+1]-R[i-1][j+1]+R[i-1][j-1]-R[i+1][j])/(4*dr*dz*(1-2*SIGMA));

				W[i][j] = -((dz^2)*(dr^2)*(1-2*SIGMA)) / ((dz^2)*(2-4*SIGMA)+(dr^2)*(4-4*SIGMA))*( H-(W[i+1][j]+W[i-1][j])/(dr^2)
					-(W[i][j+1]-W[i][j-1])/((dz^2)*(1-2*SIGMA)) );
			}
		}

		// i = 0 (left bound w/o corners)
		for (int j = 1; j < N; j++){
			R[0][j] = 0;
			W[0][j] = .5*( (2*W[1][j])/(dr^2) + (4*(1-SIGMA))/(1-2*SIGMA)*(W[0][1])/((dz^2))
						/ ((1-SIGMA)/((1-2*SIGMA)*(dz^2)) + 1/(dr^2)));
		}

		// j = N (top bound w/o corners)
		for (int i = 1; i < M; i++){
			R[i][N] = (4*R[i][N-1]-R[i][N-2]+(dz)/(dr)*(W[i-1][N]-W[i+1][N]))/3;
			W[i][N] = (4*W[i][N-1]-W[i][N-2])/(3*E*(1-SIGMA))-(2*(dz)*SIGMA)/(3*(1-SIGMA))*((R[i+1][N]-R[i-1][N])/(2*dr)+R[i][N]/r(i,N))
						+ ((P(i,N))*(1+SIGMA)*(1-2*SIGMA)*(2*dz))/(3*E*(1-SIGMA));
		}
		count++;
	}

	char a;
	std::cin >> a;
	return 0;
}
