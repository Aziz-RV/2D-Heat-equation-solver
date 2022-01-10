In order to run the code go the "build" folder, compile the src code using cmake and then the following executables are prepared:

- heat_equation_app: This is the executable for the main solver. It asks the user for the domain values (only rectangular domains), meshsize and the boundary input.
	This generates two text files (mesh.csv and numerical_solution.csv) which are stored in build/results/heat_equation_app_solutions and these can be plotted using the Matlab file "plottingfile.m".
	
- TEST_dirichlet: Runs the convergence test for the Dirichlet BC on a square grid.The results are the mesh-sizes and the corresponding errors which are displayed on the terminal.

- TEST_Neumann: Runs the convergence test for the Neumann-Dirichlett (Mixed) BC on a square grid. The results are the mesh-sizes and the corresponding errors which are displayed on the terminal.

