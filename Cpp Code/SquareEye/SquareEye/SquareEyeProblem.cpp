/*
File: SquareEyeProblem.cpp
Description: primary file for solving the Square-Eye Problem.
Author: Roland Sanford - ras9841@rit.edu
		Based on MATLAB code by John Mooney (RIT)
*/

#include <iostream>
#include <ctime>

using std::cout;


/*
Main method for running the simulation.
*/
int main(){
	// clock setup
	std::clock_t start;
	double duration;
	start = std::clock();
	//



	// display elapsed time
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	cout << "Time: " << duration << "\n";
	//
}