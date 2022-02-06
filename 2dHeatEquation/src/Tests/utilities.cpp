#include "utilities.h"


/*
    Writes the meshgrid to a .csv file which can be used to plotting
    using python or matlab.

    @params X: matrix with x-coordinates of all grid points
    @params X: matrix with y-coordinates of all grid points
    @params filename: name of file in which meshgrid is written

*/
void write_meshgrid_to_csv(Eigen::MatrixXd& X, Eigen::MatrixXd& Y,
                  std::string filename)
{

    std::ofstream csv_file(filename);

    for (int i = 0; i < X.rows();i++)
    {
        for (int j = 0; j < X.cols(); j++)
        {
            csv_file << X(i,j) << ','<< Y(i,j)<< std::endl;
        }
        
    } 
    csv_file.close();

}


/*
    Writes the result vector to a .csv file

    @params result: the result vector which is written to .csv file.
    @params filename: name of the file in which the result vector is written.
*/
void write_result_vector_to_csv(Eigen::VectorXd& result, 
                  std::string filename)
{

    std::ofstream csv_file(filename);
    for (int i = 0; i < result.size() ; i++)
    {
        csv_file << result(i) << std::endl;
    }
    csv_file.close();
}

/*
    Reads the real solution from a .csv file and writes it to a vector.

    @params filename: name of file which is read.
    @params _reference: vector on which the data of the file is written to.

*/
void read_real_solution_from_csv(Eigen::VectorXd& _reference, std::string filename )
{

    std::fstream csv_file;
    csv_file.open(filename,std::ios::in); //open the file
    int i=0;
    if (csv_file.is_open())
    {   
        std::string tp;
        while( std::getline(csv_file, tp))
        { //read data from file  and put it into string.
            double val = atof(tp.c_str()); //convert string to double
            _reference(i)=val;
            
            i++;
        }
      csv_file.close(); //close the file object.
    }
}



/*
    Discretizes the domain and each grid point is numbered in a 
    lexicohraphical fashion. It is needed for contructing other functions.

    @params Mesh: matrix which contains the numbering of each grid point.
*/
void discretize_domain(Eigen::MatrixXi& Mesh)
{   

    int Ny = Mesh.rows();
    int Nx = Mesh.cols();
    int m = 0;

    for (int i = 0; i <Ny; i++)
    {
        for (int j = 0; j< Nx; j++)
        {
            Mesh(i,j) = (i+1)*m + (j+1); 
        }
        m++;
    } 
}


/*
    Constructs the X and Y coordinate matrices which are used later to 
    plot the resulting vector in a 2D domain. 

    @params hx: meshsize along x-direction
    @params hy: meshsize along y-direction
    @params X: matrix that will contain the x-coordinates of each grid point of discretized domain.
    @params Y: matrix that will contain the y-coordinates of each grid point of doscretized domain.
*/
void meshgrid(double hx, double hy, Eigen::MatrixXd& X, 
                                    Eigen::MatrixXd& Y)
{

    int Ny = X.rows();
    int Nx = X.cols();
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j< Nx; j++)
        {
            X(i,j) = j*hx;
            Y(i,j) = i*hy;
        }
    }  
}


/*
    Takes the boundary values corresponding to the Dirichlet Boundary condition and 
    assembles the RHS vector.

    @params: u_x0: value  of u(0,y)
    @params: u_xL: value  of u(L,y)
    @params: u_y0: value  of u(x,0)
    @params: u_yL: value  of u(x,L)
    @params: RHS: vector on which the boundary values will be added.
    @params: Mesh: Matrix which contains numbering of grid points.
*/
void dirichlett_boundary(double u_x0, double u_xL, double u_y0, double u_yL
                        ,Eigen::VectorXd& RHS, Eigen::MatrixXi& Mesh)
{   
    std:: cout << "Working on RHS vector" << std::endl;
    int Nx = Mesh.cols();
    int Ny = Mesh.rows();

    int x = 0; 
    int Lx = Nx-1; 
    for (int j = 0; j < Nx; j++)
    {   
        RHS(x + j*Ny) = u_x0;   
        RHS(Lx + j*Ny) = u_xL;      
    } 
    
    int y = 0;
    int Ly = Ny-1;
    for (int i = 0; i < Ny-1; i++)
    {
        RHS(y*Nx + i) = u_y0;
        RHS(Ly*Nx + i) = u_yL;
    }

}


/*
    Takes the boundary values corresponding to the Neumann-Dirichlet Boundary condition and 
    assembles the RHS vector.

    @params: u_x0: value  of u(0,y)
    @params: u_xL: value  of du/dx(L,y)*hx (user enters du/dx(L,y))
    @params: u_y0: value  of u(x,0)
    @params: u_yL: value  of u(x,L)
    @params: RHS: vector on which the boundary values will be added.
    @params: Mesh: Matrix which contains numbering of grid points.
*/
void neumann_dirichlett_boundary(double u_x0,double u_xL, double u_y0, double u_yL ,Eigen::VectorXd& RHS, Eigen::MatrixXi& Mesh)
{   
    std:: cout << "Working on RHS vector" << std::endl;
    int Nx = Mesh.cols();
    int Ny = Mesh.rows();

    int x = 0; 
    int Lx = Nx-1; 
    for (int j = 1; j < Nx-1; j++) 
    {   
        RHS(x + j*Ny) = u_x0;      
        RHS(Lx + j*Ny) = u_xL;     
    } 
    
    int y = 0;
    int Ly = Ny-1;
    for (int i = 0; i < Ny; i++)
    {
        RHS(y*Nx + i) = u_y0;
        RHS(Ly*Nx + i) = u_yL;
    }

}




/*
    Assembles the system matrix for the system with Dirichlet boundaries. 
    The system matrix A is the discretized version of the Laplacian which accounts for the
    boundary values.
    
    @params A: Sparse system matrix
    @params Mesh: Matrix which contains numbering of grid points.
*/
void matrix_assemply_dirichlet(Eigen::SparseMatrix<double> &A, Eigen::MatrixXi& Mesh)
{   
    std:: cout << "Started Matrix Assembly" << std::endl;

    int Ny = Mesh.rows();
    int Nx = Mesh.cols();
    int N = (Nx)*(Ny);

    int Lx = Nx-1;
     A.reserve(Eigen::VectorXi::Constant(A.cols(),6));
    for (int j = 0; j< Nx; j++) 
    {
        A.coeffRef(j*Ny,j*Ny) = 1;
        A.coeffRef(Lx + j*Ny, Lx + j*Ny) = 1;
    }

    int Ly = Ny-1;
    for  (int i = 0; i < Ny-1; i++)
    {
        A.coeffRef(i,i) = 1;
        A.coeffRef(Ly*Nx + i,Ly*Nx + i) = 1;
    }

     for (int i = 0; i < N; i++)
    {
        if(A.coeffRef(i,i) != 1)
        {
            A.coeffRef(i,i) = -4;
            A.coeffRef(i,i+1) = 1;
            A.coeffRef(i,i-1) = 1;
            A.coeffRef(i, i+ Ny) = 1;
            A.coeffRef(i, i - Ny) = 1;
        }
    }
}


/*
    Assembles the system matrix for the system with mixed boundaries. 
    The system matrix A is the discretized version of the Laplacian which accounts for the
    boundary values.
    
    @params A: Sparse system matrix
    @params Mesh: Matrix which contains numbering of grid points.
*/
void matrix_assemply_neumann_dirichlet(Eigen::SparseMatrix<double> &A, Eigen::MatrixXi& Mesh)
{
    std:: cout << "Started Matrix Assembly" << std::endl;

    int Ny = Mesh.rows();
    int Nx = Mesh.cols();
    int N = (Nx)*(Ny);

    int Lx = Nx-1;
     A.reserve(Eigen::VectorXi::Constant(A.cols(),6));
    for (int j = 1; j< Nx-1; j++)
    {
        A.coeffRef(j*Ny,j*Ny) = 1;
        A.coeffRef(Lx + j*Ny, Lx + j*Ny) = 1;
        A.coeffRef(Lx + j*Ny, (Lx + j*Ny)-1) = -1;
    }

    int Ly = Ny-1;
    for  (int i = 0; i < Ny; i++)
    {
        A.coeffRef(i,i) = 1;
        A.coeffRef(Ly*Nx + i,Ly*Nx + i) = 1;
    }

     for (int i = 0; i < N; i++)
    {
        if(A.coeffRef(i,i) != 1)
        {
            A.coeffRef(i,i) = -4;
            A.coeffRef(i,i+1) = 1;
            A.coeffRef(i,i-1) = 1;
            A.coeffRef(i, i+ Ny) = 1;
            A.coeffRef(i, i - Ny) = 1;
        }
    }
}


/*
    Solves the linear sparse system using the LU solver from Eigen.

    @params A: Sparse system matrix
    @params vec: the right hand side vector of linear system.
    @return result: solution vector of linear system.
*/
Eigen::VectorXd LU_method(Eigen::SparseMatrix<double> &A,Eigen::VectorXd &vec)
{

    std::cout << "Solving system now" << std::endl;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>   solver;
    solver.analyzePattern(A); 
    solver.factorize(A); 
    Eigen::VectorXd x = solver.solve(vec);
    return x;
}


/*
    Computes the error between the numerical and real solution with respect to 2-norm

    @param x_numerical: numerical solution of pde
    @param real_x: real solution of pde
    @return : 2-norm of the difference of x_numerical and real_x
*/
double error_estimation(Eigen::VectorXd &x_numerical,Eigen::VectorXd &real_x ) // Error with 2-norm
{   
    Eigen::VectorXd error_vector= real_x-x_numerical;
    return error_vector. squaredNorm(); 
}


/*
    Computes the error between the numerical and real solution with respect to infinity-norm

    @param x_numerical: numerical solution of pde
    @param real_x: real solution of pde
    @return error: inf-norm of the difference of x_numerical and real_x
*/
double error_estimation_inf_norm(Eigen::VectorXd &x_numerical,Eigen::VectorXd &real_x )
{   
    double error = 0;
    double num;
    for (int i = 0; i<real_x.size(); i++)
    {
        num = abs(x_numerical(i)- real_x(i));
        if (num > error)
        {
            error = num;
        }
    } 
    return error;
}


/*
    Solves the linear system using the Jacobi Method. The user can enter the maximum 
    number of iterations allowed. The method stops either when the error between the 
    real solution and computed solution is smaller than a tolerance value or when the
    number of interations exceeds the maximum number of interations.

    @params A: sparse system matrix
    @params b: right hand side vector of linear system.
    @return x: solution of linear system Ax = b.
*/
Eigen::VectorXd jacobi_method(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    std::cout<< "max iterations :";
    int max_iter=0;
    std::cin >>max_iter;
    int l = b.size();
    Eigen::VectorXd x(l);
    x = Eigen::VectorXd::Zero(l);

    Eigen::VectorXd y(l);
    y = Eigen::VectorXd::Zero(l);

    double tol = 0.000001;

    double ax = 0;
    Eigen::VectorXd bb(l);
    bb = A*x;
    int k = 0;

    while(error_estimation_inf_norm(bb,b)> tol && k < max_iter)
    {
        for (int i = 0; i<l ;i++)
        {
            ax = 0;
            for (int j = 0; j < l ; j++)
            {   
                if(i != j)
                {
                    ax = ax + A.coeffRef(i,j)*x(j);
                }
            }
            y(i) = 1/A.coeffRef(i,i)*(b(i) - ax);
        }
        k = k+1;
        x = y;
        bb = A*x;
    }
    std::cout << "Total Iterates: " << k << std::endl;
    // std::cout << "Result is: " << x << std::endl;
    return x;
}


/*
    Solves the linear system using the Gauss-Seidel Method. The user can enter the maximum 
    number of iterations allowed. The method stops either when the error between the 
    real solution and computed solution is smaller than a tolerance value or when the
    number of interations exceeds the maximum number of interations.

    @params A: sparse system matrix
    @params b: right hand side vector of linear system.
    @return x: solution of linear system Ax = b.
*/
Eigen::VectorXd gauss_seidel(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    std::cout<< "max iterations :";
    int max_iter=0;
    std::cin >>max_iter;
    int l = b.size();
    Eigen::VectorXd x(l);
    x = Eigen::VectorXd::Zero(l);

    Eigen::VectorXd y(l);
    y = Eigen::VectorXd::Zero(l);

    double tol = 0.000001;

    double ax = 0;
    double ay = 0;

    Eigen::VectorXd bb(l);
    bb = A*x;
    int k = 0;

    while(error_estimation_inf_norm(bb,b)> tol && k < max_iter)
    {
        for (int i = 0; i<l ;i++)
        {
            ax = 0;
            ay = 0;
            for (int j = 0; j < i ; j++)
            {   
                ax = ax + A.coeffRef(i,j)*y(j);
            }
            for (int j = i+1; j<l ;j++)
            {
                ay = ay + A.coeffRef(i,j)*x(j);
            }
            y(i) = 1/A.coeffRef(i,i)*(b(i) - ax - ay);
        }
        k = k+1;
        x = y;
        bb = A*x;
    }
    std::cout << "Total Iterates: " << k << std::endl;
    // std::cout << "Result is: " << x << std::endl;
    return x;
}


