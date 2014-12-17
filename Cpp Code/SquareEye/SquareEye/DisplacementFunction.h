/*
File: DisplacementFunction.h
Description: header DisplacementFunction.cpp
Author: Roland Sanford - ras9841@rit.edu
*/

#ifndef DisplacementFunction_H_ 
#define DisplacementFunction_H_

class DisplacementFunction
{
	int x_val, z_val;
	char f_side[15];
	int f;

public:
	DisplacementFunction(char side[], int x, int z); 
	int getX();
	int getZ();
	int getF();
};

#endif