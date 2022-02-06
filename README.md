# PDE Solver
The project is about writing a solver for steady state solution of 2D heat equation. Since the steady state solution is time independent, the pde reduces to Poisson equation with right hand side as zero (Laplace equation) and that is what we solve in this project. We make use of the finite difference method for solving this. We plot our solution using Matlab (plotting file). We also perform convergence test for both cases (Dirichlett and Neumann-Dirichlett) to provide a proof of concept.

Already made and working:
- Solver for rectangular domain for cartesien coordinates with Dirichlett BC
- Implement a convergence test for the method with an example
- Implement Solver for Neumann-Dirichlett BC (Mixed)
- Add new linear solvers (Gauss Seidel & Jacobi)
- Organize the solver into classes
- Optimize code
- Allow negative domain values.



