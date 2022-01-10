
Here we have all the source codes.

- heat_equation.cpp: This is the main solver. This file asks the user for the domain (only rectangular domains) on which the steady statae heat equation is to be solved. Moreover, it asks the user for 			the mesh-sizes along the x and the y direction and it also asks for the boundary values. This solver works for both Dirichlet boundaries and Mixed (3 sides Dirchlet and 1 side 				Neumann) boundaries. As of now only constant boundary values are accepted. Additionally this solver allows the user to choose the linear-solver that is used to solve the linear 			system obtain upon discretizing the PDE. The three linear solvers are: 1) LU Solver (From Eigen), 2) Jacobi Solver, 3) Gauss Seidel Solver. The results of this file are two csv 			files, one of which contains the mesh (2 matrices X and Y) and the second contains the solution vector.
			** Dependencies: 1) pde_solver.h; 2) dirichlet_solver.h; 3) neumann_dirichlet_solver.h; 4) Linear_solver.h; 5) LU_solver_eigen.h; 6) Jacobi_solver.h; 7) Gauss_Seidel_solver.h; 
			8) Eigen**
			
			
//////////////////////////////////////////////////////////////////////////////////////////////

In src/solvers/ we have all the dependencies of heat_equation.cpp

This includes: 

- pde_solver.h that disceretizes the pde, creates a mesh and solver the linear system.

- neumann_dirichlet_solver.h: prepares the system matrix and the rhs vector for the mixed boundary system.

- dirichlet_solver.h: prepares the system matrix and the rhs vector for the dirichlet boundary system.

In /solvers/linear_solvers/ we have the linear iterative solvers that solver the sparse linear system

- Linear_solver.h: Is the base class of solvers. This takes the sparse matrix (system matrix), the rhs vector and the result vector (which is computed by the linear_solvers)

- LU_solver_eign.h: Is a derived class of Linear_solver and uses the LU solver from the library Eigen to solve the discretized linear system.

- jacobi_solver.h: Is a derived class of Linear_solver and uses the Jacobi Method to solve the discretized linear system.

- gauss_seidel_solver.h: Is a derived class of Linear_solver and uses the Gauss Seidel to solve the discretized linear system.


/////////////////////////////////////////////////////////////////////////////////////////////

In src/Tests/ we have
- utilities.cpp: contains the functions that are used by the convergence tests "test_dirichlet.cpp" and "test_neumann.cpp".
		 ** Dependencies: 1) Eigen **


- test_Dirichlet.cpp: This file runs the convergence test of the method for a specific test case with Dirichlet boundaries on a square domain. This file generates the solution vectors for 4 different 			mesh-sizes and compares the numerical solution to the reference solution files in the \build folder. And the end the file prints the error corresponding to each of the mesh sizes. 
			Please note that this code contains an extra function called 'dirichlett_boundary_function' that is not used in the main solver (heat_equation.cpp) as for this case
			one of the boundary values is actually a function of x-coordinate. This has been done purposely such that the analytical solution has a form that not an eigenfunction-expansion as 				it would be in case of constant boundary values.
			** Dependencies: 1) utilities.h ; 2) math.h **

- test_neumann.cpp: This file runs the convergence test of the method for a specific test case with neumann boundaries on a square domain. This file generates the solution vectors for 4 different 				mesh-sizes and compares the numerical solution to the reference solution files in the \build folder. And the end the file prints the error corresponding to each of the mesh sizes. 
			Please note that this code contains an extra function called 'neumann_boundary_function' that is not used in the main solver (heat_equation.cpp) as for this case
			one of the boundary values is actually a function of x-coordinate. This has been done purposely such that the analytical solution has a form that not an eigenfunction-expansion as 				it would be in case of constant boundary values.
			** Dependencies: 1) utilities.h ; 2) math.h **			
			
Details of the convergence tests are in the README file in src/Tests/.

////////////////////////////////////////////////////////////////////////////////////////////








