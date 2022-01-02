#include "utilities.h" 
#include <math.h>
// we can add all the includes but they are already in utilities.h



Eigen::VectorXd neumann_boundary_function(int N,double h, double u_xL, double u_y0, 
                                            Eigen::MatrixXd& X,Eigen::MatrixXd& Y)
{
    /*
        This function is only necessary for this test case where a boundary prescription 
        is a function and not a value.
    */
    Eigen::VectorXd RHS(N*N);
    RHS= Eigen::VectorXd::Zero(N*N);//we don't know we are intializing RHS to 0 but this way everything is working
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
        u_yL = cos(pi*num)*std::sinh(pi); // since u(x,L) = sin(pi*x)*sinh(pi) the values of x-coordinate are needed
                                        // and this is just any row of the X-matrix, values of Y are not needed.  
        RHS(y*N + i) = u_y0;

        RHS(Ly*N + i) = u_yL;
    }
    return RHS;
}

void pde_neumann_solver_test(double h, 
                                    std::string ref_filename,std::string solution_filename)
{

    int x0 = 0;
    int xL = 1;
    double hx = h; // mesh size
    int Nx = (xL-x0)/hx +1; // number of columns
    

    int y0 = 0;
    int yL = 1;
    double hy = h; // mesh size
    int Ny = (yL-y0)/hy +1; // No. of rows
    // Meshing
    Eigen::MatrixXi Mesh(Ny,Nx); 
    discretize_domain(Mesh);

    Eigen::MatrixXd X(Nx,Ny); 
    Eigen::MatrixXd Y(Nx,Ny);
    meshgrid(hx,hy,X,Y); 

    write_meshgrid_to_csv(X,Y);

    Eigen::VectorXd RHS(Nx*Ny);
    
    double u_xL = 0; // du/dx (L,y) = 0
    double u_y0 = 0; // u(x,0) = 0
    
    RHS=neumann_boundary_function(Nx,h, u_xL ,u_y0, X, Y);

    Eigen::SparseMatrix<double> A(Nx*Ny, Nx*Ny);
    matrix_assemply_neumann_dirichlet(A, Mesh);
    

    Eigen::VectorXd U(Nx*Ny);
    U=LU_method(A,RHS);
  
    
    // Save result to csv
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
        We are going to make a convergence test with these values of step size: (0.2,0.1,0.05,0.025)
        and check if it converges and if the error is divided by 4 for every halving.   
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
