#include "utilities.h"


void write_meshgrid_to_csv(Eigen::MatrixXd& X, Eigen::MatrixXd& Y,
                  std::string filename)
{
    /*
        Writes the meshgrid to a file which can be used to plotting
        using python or matlab.
    */
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


void write_result_vector_to_csv(Eigen::VectorXd& result, 
                  std::string filename)
{
    /*
        Writes the solution vector to a csv file.
    */
    std::ofstream csv_file(filename);
    for (int i = 0; i < result.size() ; i++)
    {
        csv_file << result(i) << std::endl;
    }
    csv_file.close();
}


void read_real_solution_from_csv(Eigen::VectorXd& _reference, std::string filename )
{
    /*
        Reads the real solution from a csv file 
    */
    std::fstream csv_file;
    csv_file.open(filename,std::ios::in); //open the file
    int i=0;
    if (csv_file.is_open())
    {   //checking whether the file is open
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



void discretize_domain(Eigen::MatrixXi& Mesh)
{   
    /*
        Discretizes the domain and each grid point is numbered in a 
        lexicohraphical fashion. It is needed for contructing other functions.
    */
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



void meshgrid(double hx, double hy, Eigen::MatrixXd& X, 
                                    Eigen::MatrixXd& Y)
{
    /*
        Constructs the X and Y coordinate matrices which are used later to 
        plot the resulting vector in a 2D domain. 
    */
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



void dirichlett_boundary(double u_x0, double u_xL, double u_y0, double u_yL
                        ,Eigen::VectorXd& RHS, Eigen::MatrixXi& Mesh)
{   
    std:: cout << "Working on RHS vector" << std::endl;
    int Nx = Mesh.cols();
    int Ny = Mesh.rows();

    int x = 0; // First Column
    int Lx = Nx-1; // Last column
    for (int j = 0; j < Nx; j++)
    {   
        RHS(x + j*Ny) = u_x0;       // Makes First column 0 // u(0,y) =  0
        //std::cout<< x + j*Ny << ':' << u_x0 << std::endl;
        RHS(Lx + j*Ny) = u_xL;      // Makes Last Column 0 // u(L,y) = 0
        //std::cout<< Lx + j*Ny << ':' << u_xL << std::endl;
    } 
    
    int y = 0;
    int Ly = Ny-1;
    for (int i = 0; i < Ny-1; i++)
    {
        RHS(y*Nx + i) = u_y0;
        //std::cout<< y*Nx + i << ':' << u_y0 << std::endl;
        RHS(Ly*Nx + i) = u_yL;
        //std::cout<< Ly*Nx + i << ':' << u_yL << std::endl;
    }

}




void neumann_dirichlett_boundary(double u_xL,double u_x0, double u_y0, double u_yL ,Eigen::VectorXd& RHS, Eigen::MatrixXi& Mesh)
{   
    std:: cout << "Working on RHS vector" << std::endl;
    int Nx = Mesh.cols();
    int Ny = Mesh.rows();

    int x = 0; // First Column
    int Lx = Nx-1; // Last column
    for (int j = 1; j < Nx-1; j++) // changed things here
    {   
        RHS(x + j*Ny) = u_x0;       // Makes First column 0 // u(0,y) =  0
        //std::cout<< x + j*Ny << ':' << u_x0 << std::endl;
        RHS(Lx + j*Ny) = u_xL;      // Makes Last Column L // (du/dx)(L,y) = g(j)
        //std::cout<< Lx + j*Ny << ':' << u_xL << std::endl;
    } 
    
    int y = 0;
    int Ly = Ny-1;
    for (int i = 0; i < Ny; i++)
    {
        RHS(y*Nx + i) = u_y0;
        //std::cout<< y*Nx + i << ':' << u_y0 << std::endl;
        RHS(Ly*Nx + i) = u_yL;
        //std::cout<< Ly*Nx + i << ':' << u_yL << std::endl;
    }

}

//not used
void neumann_boundary(double du_xL,double du_x0, double du_y0, double du_yL ,Eigen::VectorXd& RHS, Eigen::MatrixXi& Mesh)
{   
    std:: cout << "Working on RHS vector" << std::endl;
    int Nx = Mesh.cols();
    int Ny = Mesh.rows();

    int x = 0; // First Column
    int Lx = Nx-1; // Last column
    for (int j = 1; j < Nx-1; j++) // changed things here
    {   
        RHS(x + j*Ny) = du_x0;       // Makes First column 0 // du/dx(0,y) =  g(1)
        //std::cout<< x + j*Ny << ':' << u_x0 << std::endl;
        RHS(Lx + j*Ny) = du_xL;      // Makes Last Column L // (du/dx)(L,y) = g(2)
        //std::cout<< Lx + j*Ny << ':' << u_xL << std::endl;
    } 
    
    int y = 0;
    int Ly = Ny-1;
    for (int i = 0; i < Ny; i++)
    {
        RHS(y*Nx + i) = du_y0;// du/dy(x,0) =  g(3)
        //std::cout<< y*Nx + i << ':' << u_y0 << std::endl;
        RHS(Ly*Nx + i) = du_yL;//du/dy(x,0) =  g(3)
        //std::cout<< Ly*Nx + i << ':' << u_yL << std::endl;
    }

}


void matrix_assemply_dirichlet(Eigen::SparseMatrix<double> &A, Eigen::MatrixXi& Mesh)
{   
    std:: cout << "Started Matrix Assembly" << std::endl;
    /*
        This function constructs system matrix for the discretized system.
        This system matrix takes depends on the boundary conditions and is the discretized form
        of the laplacian operator.
    */
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

void matrix_assemply_neumann_dirichlet(Eigen::SparseMatrix<double> &A, Eigen::MatrixXi& Mesh)
{
    std:: cout << "Started Matrix Assembly" << std::endl;
    /*
        This function constructs system matrix for the discretized system.
        This system matrix takes depends on the boundary conditions and is the discretized form
        of the laplacian operator.
    */
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
//not used 
void matrix_assemply_neumann(Eigen::SparseMatrix<double> &A, Eigen::MatrixXi& Mesh)
{
    std:: cout << "Started Matrix Assembly" << std::endl;
    /*
        This function constructs system matrix for the discretized system.
        This system matrix takes depends on the boundary conditions and is the discretized form
        of the laplacian operator.
    */
    int Ny = Mesh.rows();
    int Nx = Mesh.cols();
    int N = (Nx)*(Ny);

    int Lx = Nx-1;
     A.reserve(Eigen::VectorXi::Constant(A.cols(),6));
    for (int j = 1; j< Nx-1; j++)
    {
        A.coeffRef(j*Ny,j*Ny) = 1;
        A.coeffRef(j*Ny,j*Ny+1) = -1;
        A.coeffRef(Lx + j*Ny, Lx + j*Ny) = 1;
        A.coeffRef(Lx + j*Ny, (Lx + j*Ny)-1) = -1;
    }

    int Ly = Ny-1;
    for  (int i = 0; i < Ny; i++)
    {
        A.coeffRef(i,i) = 1;
        A.coeffRef(i,i+Nx) = -1;

        A.coeffRef(Ly*Nx + i,Ly*Nx + i) = 1;
        A.coeffRef(Ly*Nx + i,Ly*Nx + i-Nx) = -1;
        

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


Eigen::VectorXd LU_method(Eigen::SparseMatrix<double> &A,Eigen::VectorXd &vec)
{
    /*
        Solves the linear system to obtain the solution vector 
    */
    std::cout << "Solving system now" << std::endl;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>   solver;
    solver.analyzePattern(A); 
    solver.factorize(A); 
    Eigen::VectorXd x = solver.solve(vec);
    return x;
}

double error_estimation(Eigen::VectorXd &x_numerical,Eigen::VectorXd &real_x ) // Error with 2-norm
{   
    /*
        Computes the 2-norm of the error in solution vector and the real solution 
    */
    Eigen::VectorXd error_vector= real_x-x_numerical;
    return error_vector. squaredNorm(); 
}

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
