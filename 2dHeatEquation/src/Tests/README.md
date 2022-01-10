
		      
Details of the Dirichlet Convergence test

The convergence test is meant to check the validity of our implementation. 
The finite difference method solves the steady state heat equation on a 2D square domain.
The test case is on domain [0,1]^2 with meshsize being equal along both x and y direction.
Results of convergence test are with respect to inf_norm.
The boundary conditions are as follows:

u(x,0) = 0,

u(0,y) = 0,

u(y,1) = 0,

u(x,1) = sin(pi * x)*sinh(pi)

Analytical solution: u(x,y) = sin(pi * x)*sinh(pi * y)


//////////////////////////////////////////////////////////////////////////////


Details of the Neumann-Dirichlet Convergence Test

The convergence test is meant to check the validity of our implementation. 
The finite difference method solves the steady state heat equation on a 2D square domain.
The test case is on domain [0,1]^2 with meshsize being equal along both x and y direction.
Results of convergence test are with respect to inf_norm.
The boundary conditions are as follows:


u(x,1) = cos(pi * x)* sinh(pi)


u(0,y) = sinh(pi * y)


u(x,0) = 0


du(x,y)/dx = 0 (at x = 1)


Analytical Solution: u(x,y) = cos(pi * x)* sinh(pi* y)

///////////////////////////////////////////////////////////////////////////////
