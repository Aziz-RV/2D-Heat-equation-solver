
Here we have all the source codes.
- utilities.cpp: contains the functions that are used by both test_dirichlett.cpp and heat_equation.cpp

- heat_equation.cpp: This is the main solver.
		      This file asks the user for the domain on which the heat equation is to be solved. Moreover, it asks the user for the mesh-sizes along the x and the y direction
		      and it also asks for the boundary values. 
		      Note that this code only works for constant Dirichlett boundary values. 
		      The results of this file are two csv files, one of which contains the mesh (2 matrices X and Y) and the second contains the solution vector.

- test_dirichlett.cpp: This file runs the convergence test of the method for a specific test case. This file generates the solution vectors for 4 different mesh
			sizes and compares the numerical solution to the reference solution files in the \build folder. And the end the file prints the error corresponding
			to each of the mesh sizes. 
			Please note that this code contains an extra function called 'dirichlett_boundary_function' that is not used in the main solver as for this case
			one of the boundary values is actually a function of x-coordinate. This has been done purposely such that the analytical solution has a form that
			not an eigenfunction-expansion as it would be in case of constant boundary values.
			
			
Details of the Convergence test

The convergence test is meant to check the validity of our implementation. 
The finite difference method solves the steady state heat equation on a 2D square domain.
The test case is on domain [0,1]^2 with meshsize being equal along both x and y direction.
The boundary conditions are as follows:

u(x,0) = 0,

u(0,y) = 0,

u(y,1) = 0,

u(x,1) = sin(pi * x)*sinh(pi)

the analytical solution is u(x,y) = sin(pi * x)*sinh(pi * y)




