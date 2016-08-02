<a name="clp">ContactLensProject</a>
=====================================
Maintained by Roland Sanford (<ras9841@rit.edu>).
  
##Table of Contents

1. [Setup](#1)
  1. [Development Conditions](#1.1)
  2. [Compiling and Running](#1.2)
2. [Implementation](#2)
  1. [Gauss-Seidel method](#2.1)
  2. [Finite-difference method](#2.2)
  3. [Computational Grid](#2.3)
  4. [General pseudocode](#2.4)

##<a name="1"></a>1. Setup [[top](#clp)]
###<a name="1.1"></a>1.1 Development Conditions
The code for this project was written in C++ using the gcc (GCC) 4.9.2 (Red Hat 4.9.2-6) GNU compiler.
This repository is under active development, so the current repo head may not compile.

###<a name="1.2"></a>1.2 Compiling and Running
The information below details how to compile and run the project's code for various operating systems.
You can either clone or download the repository (and `build`, if you are running Linux).

####Linux
Open a terminal and navigate to the main project directory. Type 
```{r, engine='bash'}
$ g++ -o <executable name> c_eye.cpp helper.cpp
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
$ clang++ -std=c++11 -stdlib=libc++ c_eye.cpp helper.cpp
$ ./a.out
```
to compile then run the program.

###Windows
This requires `Visual C++ .NET`. Run Visual C++ 2010 Express Command Prompt. Type
```{r, engine='bash'}
> cl /EHsc c_eye.cpp helper.cpp
> c_eye.exe
```
to compile then run the program.

##<a name="2"></a>2. Implementation [[top](#clp)]
###<a name="2.1"></a>2.1 Gauss-Seidel method
We use the Gauss-Seidel method (GSM) is our general technique for solving the cylindrical 
eye problem. In this section, we will provide a brief introduction to the method.  
  
GSM is used to solve linear systems of equations in the form `Ax=b` one equation 
at a time. During each iteration, the new values are calculated based on the values 
from the previous iteration. An initial guess serves as the "old" values for the 
first run through. This process continues continues until a predefined convergence 
condition is satisfied. GSM convergence requires the matrix `A` to be strictly 
diagonally dominant, or symmetric and positive definite. 
    
###<a name="2.2"></a>2.2 Finite-difference method [[top](#clp)]
We discretized the governing partial differntial equations (PDEs) and boundary conditions 
using finite-difference method (FDM). FDM uses the values of a point and it's neighbors 
to approximate derivatives at that point. Combinations of forward, backward, and central 
finite difference equations were used thoughout the computational grid. Currently,
a global second-order accuracy is held throughout the simulation.
  
###<a name="2.3"></a>2.3 Computational grid [[top](#clp)]
We reprensent the eye as a radially-symmetric cylinder. Consequently, we only need 
sole the domain from the center of the eye to an edge to solve the system. The eye 
is mapped to an `MxN` grid. The variable `i` is used to represent rows and height, 
such that `0<=i<=M`, and `j` to represent columns and width, such that `0<=j<=N`. 
As a whole, the grid can be viewed as follows:
```
[M][0]  |---------------| [M][N]
 ^      |---------------|
 i      |---------------|
 v      |---------------|
[0][0]  |---------------| [0][N]
        <------ j ------>
```
Note: `[]` denotes array indexing in C (zero-based). 
  
###<a name="2.4"></a>2.4 General pseudocode [[top](#clp)]
Here, we present a high-level overview of the simulation. 
```
Initialize the displacement functions R and W to be M+1xN+1 arrays.
Make the initial guess: R=0, W=0.
While the system has not converged to solution:
    Every k iterations update the suction pressure.
	Solve for R and W at each point on the computational grid.
```

