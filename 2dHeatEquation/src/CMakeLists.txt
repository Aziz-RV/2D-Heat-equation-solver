cmake_minimum_required (VERSION 3.0)
project (heat_equation_app VERSION 0.2.0 LANGUAGES CXX)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#heat_equation is for the user use
add_executable (heat_equation_app heat_equation.cpp
 solvers/pde_solver.cpp solvers/linear_solvers/Linear_solver.cpp solvers/linear_solvers/gauss_seidel_solver.cpp solvers/linear_solvers/jacobi_solver.cpp solvers/linear_solvers/LU_solver_eigen.cpp 
 solvers/dirichlet_solver.cpp solvers/neumanndirichlet_solver.cpp)
target_link_libraries (heat_equation_app Eigen3::Eigen)
# We also want to build a separate executable with our tests.
# We have moved the test of Dirichilet boundaries into Test_Dirichlet.cpp
add_executable (TEST_dirichlet Tests/test_Dirichlet.cpp Tests/utilities.cpp)
# Add a test case, which simply runs some tests with different configurations
target_link_libraries (TEST_dirichlet Eigen3::Eigen)
add_executable (TEST_Neumann Tests/test_neumann.cpp Tests/utilities.cpp)
# Add a test case, which simply runs some tests with different configurations
target_link_libraries (TEST_Neumann Eigen3::Eigen)
