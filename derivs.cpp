//
// File:	derivs.cpp
// 
// Implements forwards, backwards, and centered finite difference 
// equations to 4th order. Also contains derivative operators for 
// the mapped computational coordinates (r_hat and z_hat).
//
// @author Roland Sanford <ras9841.rit.edu>

#include "derivs.h"

// c_diff1 computes a 4th order centered finite difference of the 
// first derivative of the function f about f[r][c].
//
// Preconditions:
//		f -> function to be operated on having a spread of 2 indecies.
//		r -> row index of the function at the derivative location.
//		c -> column index of the function at the derivative location.
//		dir -> direction of the derivative (x or y)
//		ds -> distance between two points in the dir direction.
// Postconditions:
//		The derivative at the point f[r][c] has been calculated to 
//		4th order accuracy.
double c_diff1(double **f, size_t r, size_t c, char dir, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal derivative
	{
		deriv = (-f[r+2][c] + 8*f[r+1][c] - 8*f[r-1][c] + f[r-2][c])
					/ (12*ds);
	}
	else // Verticle derivative
	{
		deriv = (-f[r][c+2] + 8*f[r][c+1] - 8*f[r][c-1] + f[r][c-2])
					/ (12*ds);
	}
	return deriv;
}

// f_diff1 computes a 4th order forward finite difference of the 
// first derivative of the function f about f[r][c].
//
// Preconditions:
//		f -> function to be operated on at left bound.
//		r -> row index of the function at the derivative location.
//		c -> column index of the function at the derivative location.
//		dir -> direction of the derivative (x or y)
//		ds -> distance between two points in the dir direction.
// Postconditions:
//		The derivative at the point f[r][c] has been calculated to 
//		4th order accuracy.
double f_diff1(double **f, size_t r, size_t c, char dir, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal derivative
	{
		deriv = (f[r+4][c] - 6*f[r+3][c] + 18*f[r+2][c] - 10*f[r+1][c]
				- 3*f[r][c])/(12*ds);
	}
	else // Verticle derivative
	{
		deriv = (f[r][c+4] - 6*f[r][c+3] + 18*f[r][c+2] - 10*f[r][c+1]
				- 3*f[r][c])/(12*ds);
	}
	return deriv;
}

// b_diff1 computes a 4th order backward finite difference of the 
// first derivative of the function f about f[r][c].
//
// Preconditions:
//		f -> function to be operated on at left bound.
//		r -> row index of the function at the derivative location.
//		c -> column index of the function at the derivative location.
//		dir -> direction of the derivative (x or y)
//		ds -> distance between two points in the dir direction.
// Postconditions:
//		The derivative at the point f[r][c] has been calculated to 
//		4th order accuracy.
double b_diff1(double **f, size_t r, size_t c, char dir, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal derivative
	{
		deriv = (3*f[r][c] + 10*f[r-1][c] - 18*f[r-2][c] + 6*f[r-3][c]
				- f[r-4][c])/(12*ds);
	}
	else // Verticle derivative
	{
		deriv = (3*f[r][c] + 10*f[r][c-1] - 18*f[r][c-2] + 6*f[r][c-3]
				- f[r][c-4])/(12*ds);	
	}
	return deriv;
}

