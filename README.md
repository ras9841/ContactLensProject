<a name="clp">ContactLensProject</a>
=====================================
Maintained by Roland Sanford (<ras9841@rit.edu>).
<a name="top"></a>
##Table of Contents

1. [Setup](#1)
  1. [Development Conditions](#1.1)
  2. [Compiling and Running](#1.2)
2. [Implementation](#2)
  1. [Gauss-Seidel method](#2.1)
  2. [Finite difference method](#2.2)


##<a name="1"></a>1. Setup [[top](#clp)]
###<a name="1.1"></a>1.1 Development Conditions
The code for `c_eye.cpp` was written in C++ using the gcc (GCC) 4.9.2 (Red Hat 4.9.2-6) GNU compiler.
 
###<a name="1.2"></a>1.2 Compiling and Running
The information below details how to compile and run the `c_eye.cpp` file for various operating systems.
You can either clone the repository or download `c_eye.cpp` (and `build`, if you are running Linux).

####Linux
Open a terminal and navigate to the main project directory. Type 
```{r, engine='bash'}
$ g++ -o <executable name> c_eye.cpp
```
to manually compile the program. To run, enter
```{r, engine='bash'}
$ ./<executable name>
```

Alternatively, your can run `build` to produce an executable calld `clp`. This is done in the following way:
```{r, engine='bash'}
$ sh build 
$ ./clp
```

###Mac
Open a terminal and navigate to the main project directory. This requires XCode. Type 
```{r, engine='bash'}
$ clang++ -std=c++11 -stdlib=libc++ c_eye.cpp
$ ./a.out
```
to compile then run the program.

###Windows
This requires `Visual C++ .NET`. Run Visual C++ 2010 Express Command Prompt. Type
```{r, engine='bash'}
> cl /EHsc c_eye.cpp
> c_eye.exe
```
to compile then run the program.

##<a name="2"></a>2. Implementation [[top](#clp)]
###<a name="2.1"></a>2.1 Gauss-Seidel method
The general technique we used to solve our cylindrical eye problem was the Gauss-Seidel method (GSM). Here, we will provide a brief introduction to the method.  
  
GSM is used to solve linear systems of equations in the form `Ax=b` one equation at a time. During each iteration, the new values are calculated from the old values. An initial guess is initialized to serve as the "old" values for the first run through. This process continues continues until a predefined convergence condition is satisfied. GSM requires the matrix `A` to be strictly diagonally dominant or symmetric positive. 
  
###<a name="2.2"></a>2.2 Finite difference method

