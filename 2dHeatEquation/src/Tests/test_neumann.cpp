#include "utilities.h" 
#include <math.h>

/*
        This function is only necessary for this test case where a boundary prescription 
        is a function and not a value.

        This function assembles the RHS vector when the boundary value has a specific functional value.

        @params N: number of elements in right hand side vector
        @params h: meshsize
        @params u_xL: neumann boundary value ((du/dx)*h = u_xL)
        @params u_y0: dirichlet boundary valye u(x,0)
        @params X: matrix with x-coordinates of grid points
        @params Y: matrix with y-coordinates of grid points
        @return RHS: the right hand side vector of discretized linear system which accounts for boundary values
*/
Eigen::VectorXd neumann_boundary_function(int N,double h, double u_xL, double u_y0, 
                                            Eigen::MatrixXd& X,Eigen::MatrixXd& Y)
{

    Eigen::VectorXd RHS(N*N);
    RHS= Eigen::VectorXd::Zero(N*N);
    double num;
    double u_yL;
    double u_x0;
    const double pi = 3.141592653589;
    int x = 0; 
    int Lx = N-1; 
    for (int j = 1; j < N-1; j++)
    {   
        num = Y(j,0);
        u_x0 = std::sinh(pi*num);
        RHS(x + j*N) = u_x0; 
        RHS(Lx + j*N) = u_xL*h;   //neumann value
    } 
    
    int y = 0;
    int Ly = N-1;

    for (int i = 0; i < N; i++)
    {

        num = X(0,i);
        u_yL = cos(pi*num)*std::sinh(pi); 
                                         
        RHS(y*N + i) = u_y0;

        RHS(Ly*N + i) = u_yL;
    }
    return RHS;
}


/*
    Performs all the necessary tasks to solve the steady state heat equation. 
    - Meshing (discretizing)
    - Preparing Meshgrid (& writing it to a csv file)
    - Assembling RHS vector
    - Assembling system matrix
    - Solving linear system using LU_Solver (from Eigen) & writing it to a csv file.
    - Reading the reference files to obtain the reference solution
    - Compute inf-norm between the reference and the numerical solution.

    @params h: meshsize
    @params ref_filename: name of the file the contains the reference solution vector
    @params solution_filename": name of the file in which the solution vector will be written.
*/
void pde_neumann_solver_test(double h, 
                                    std::string ref_filename,std::string solution_filename)
{

    int x0 = 0;
    int xL = 1;
    double hx = h; // mesh size
    int Nx = (xL-x0)/hx +1; // number of columns
    

    int y0 = 0;
    int yL = 1;
    double hy = h; 
    int Ny = (yL-y0)/hy +1; 

    Eigen::MatrixXi Mesh(Ny,Nx); 
    discretize_domain(Mesh);

    Eigen::MatrixXd X(Nx,Ny); 
    Eigen::MatrixXd Y(Nx,Ny);
    meshgrid(hx,hy,X,Y); 

    write_meshgrid_to_csv(X,Y);

    Eigen::VectorXd RHS(Nx*Ny);
    
    double u_xL = 0; 
    double u_y0 = 0; 
    
    RHS=neumann_boundary_function(Nx,h, u_xL ,u_y0, X, Y);

    Eigen::SparseMatrix<double> A(Nx*Ny, Nx*Ny);
    matrix_assemply_neumann_dirichlet(A, Mesh);
    
    Eigen::VectorXd U(Nx*Ny);
    U=LU_method(A,RHS);
  
    write_result_vector_to_csv(U,solution_filename);
    Eigen::VectorXd _reference(Nx*Ny);
    read_real_solution_from_csv(_reference, ref_filename);

    double error;
    error = error_estimation_inf_norm(U,_reference); 
    std::cout << "for the step size :"<<h << " we have an error of :"<< error << std::endl;
    
}


int main()
{
    /*
        This runs a convergence test for the Dirichlet boundaries for the steady state heat equation on a square domain.
        The error should become halved when the meshsize is halved. 
        This indicates that the method has O(1) convergence as we used a first order difference method 
        for Neumann boundary.
    */
    Eigen::Vector4d h(0.2,0.1,0.05,0.025); // Mesh sizes for convergence test
    std::cout << "h , Error" << std::endl;

    pde_neumann_solver_test(h(0),"./testneumann_references_data/reference1.csv","./results/test_neumann/sol_1.csv");
    std::cout <<std::endl;
    pde_neumann_solver_test(h(1),"./testneumann_references_data/reference2.csv","./results/test_neumann/sol_2.csv");
    std::cout << std::endl;
    pde_neumann_solver_test(h(2),"./testneumann_references_data/reference3.csv","./results/test_neumann/sol_3.csv");
    std::cout <<std::endl;
    pde_neumann_solver_test(h(3),"./testneumann_references_data/reference4.csv","./results/test_neumann/sol_4.csv");
    

    return 0;
}
