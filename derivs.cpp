//
// File:	derivs.cpp
// 
// Implements forwards, backwards, and centered finite difference 
// equations to 2nd order. Also contains derivative operators for 
// the mapped computational coordinates (r_hat and z_hat).
//
// @author Roland Sanford <ras9841.rit.edu>

#include "derivs.h"

// c_diff1 computes a 4th-order centered finite difference of the 
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
//		4th-order accuracy.
double c_diff1(double **f, char dir, size_t r, size_t c, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal
	{
		deriv=(-f[r][c+2]+f[r][c-2]+8*(f[r][c+1]-f[r][c-1]))/(12*ds);
	}
	else // Verticle
	{
		deriv=(-f[r+2][c]+f[r-2][c]+8*(f[r+1][c]-f[r-1][c]))/(12*ds);
	}
	return deriv;
}

// c_diff2 computes a 4th-order centered finite difference of the 
// second derivative of the function f about f[r][c].
//
// Preconditions:
//		f -> function to be operated on having a spread of 2 indecies.
//		r -> row index of the function at the derivative location.
//		c -> column index of the function at the derivative location.
//		dir -> direction of the derivative (x or y)
//		ds -> distance between two points in the dir direction.
// Postconditions:
//		The second derivative at the point f[r][c] has been calculated 
//		to 4th-order accuracy.
double c_diff2(double **f, size_t r, size_t c, char dir, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal
	{
		deriv = (-(f[r][c]
	}
	else // Verticle
	{i
		deriv = (f[r+1][c] - 2*f[r][c] + f[r-1][c])/(ds^2);
	}
	return deriv;
}

// f_diff1 computes a 2nd-order forward finite difference of the 
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
//		2nd-order accuracy.
double f_diff1(double **f, size_t r, size_t c, char dir, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal derivative
	{
		del1 = f[r][c+1]-f[r][c];
		del2 = f[r][c+2]-f[r][c];
		del3 = f[r][c+3]-f[r][c];

		deriv = (3*del1-3*del2/2.0+del3/3.0)/(ds);
	}
	else // Verticle
	{
		del1 = f[r+1][c]-f[r][c];
		del2 = f[r+2][c]-f[r][c];
		del3 = f[r+3][c]-f[r][c];

		deriv = (3*del1-3*del2/2.0+del3/3.0)/(ds);
	}
	return deriv;
}

// f_diff2 computes a 2nd-order forward finite difference of the 
// second derivative of the function f about f[r][c].
//
// Preconditions:
//		f -> function to be operated on at left bound.
//		r -> row index of the function at the derivative location.
//		c -> column index of the function at the derivative location.
//		dir -> direction of the derivative (x or y)
//		ds -> distance between two points in the dir direction.
// Postconditions:
//		The second derivative at the point f[r][c] has been calculated to 
//		2nd-order accuracy.
double f_diff2(double **f, size_t r, size_t c, char dir, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal derivative
	{
		del1 = f[r][c+1]-f[r][c];
		del2 = f[r][c+2]-f[r][c];
		del3 = f[r][c+3]-f[r][c];

		deriv = (-5*del1+4*del2-del3)/(dx^2);
	}
	else // Verticle
	{
		del1 = f[r+1][c]-f[r][c];
		del2 = f[r+2][c]-f[r][c];
		del3 = f[r+3][c]-f[r][c];

		deriv = (-5*del1+4*del2-del3)/(dx^2);
	}
	return deriv;
}

// b_diff1 computes a 2nd-order backward finite difference of the 
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
//		2nd order accuracy.
double b_diff1(double **f, size_t r, size_t c, char dir, double ds)
{
	double deriv;
	if (dir == 'x') // Horizontal derivative
	{
		del1 = f[r][c]-f[r][c-1];
		del2 = f[r][c]-f[r][c-2];
		del3 = f[r][c]-f[r][c-3];

		deriv = (-3*del1+3*del2/2.0-del3/3)/(ds);
	}
	else // Verticle derivative
	{
	}
	return deriv;
}

