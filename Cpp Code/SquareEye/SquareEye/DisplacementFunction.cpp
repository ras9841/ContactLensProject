/*
File: DisplacementFunction.cpp
Description: blueprint for a displacement function U (x-y direction) and W (z direction)
Author: Roland Sanford - ras9841@rit.edu
		Based on MATLAB code by John Mooney (RIT)
*/

#include "DisplacementFunction.h"
#include <stdexcept>

/*
Constructor for the displacement function.
*/

DisplacementFunction::DisplacementFunction(char side[], int x, int z){
	x_val = x;
	z_val = z;
	for (int i = 0; side[i] != '\0'; i++){
		f_side[i] = side[i];
	}

	switch(f_side[0]){
	case 'l':
	case 'r':
	case 'b':
		f = 0;
		break;
	deafult:
		throw std::invalid_argument("Received invalid side.");
	}
}

int DisplacementFunction::getX(){
	return x_val;
}

int DisplacementFunction::getZ(){
	return z_val;
}

int DisplacementFunction::getF(){
	return f;
}