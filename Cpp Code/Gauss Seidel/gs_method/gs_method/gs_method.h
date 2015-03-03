//
// File: gs_method.h
//
// Header file for the Gauss-Seidel method. Allows for the sharing
// of equation variables and of functions.
//
// @author Roland Sanford <ras9841@rit.edu>
// 
// // // // // // // // // // // // // // // // // // // // // // // // //

#ifndef _GS_METHOD_H
#define _GS_METHOD_H

#include <vector>

// 
// Globals
//

extern const int M,N;
extern double E, SIGMA;

//
// Shared Functions
//
void gs_method(int i, int j, std::vector<std::vector<double>>& R, std::vector<std::vector<double>>& W);

#endif