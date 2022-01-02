#include "utilities.h" 
// we can add all the includes but they are already in utilities.h


int main()
{
    std::cout << "Enter the domain:" << std::endl;
    std::cout << "Enter the initial and final x-value:" << std::endl;
    double x0;
    double xL;
    std::cout << "Enter the initial of x:" ;
    std::cin >> x0;
    std::cout << "and the final value of x:" ;
    std::cin >> xL;

   
    double y0;
    double yL;
    std::cout << "Enter the initial of y:" ;
    std::cin >> y0;
    std::cout << "and the final value of y:" ;
    std::cin >> yL;


    double hx; // mesh size
    double hy; // mesh size
    std::cout << "Enter the mesh size in x and y:" << std::endl;
    std::cout<<"hx= ";std::cin >> hx;
    std::cout<<"hy= ";std::cin >> hy;
    int Nx = (xL-x0)/hx +1; // No. of columns (must be int)
    int Ny = (yL-y0)/hy +1; // No. of rows


    // Meshing
    Eigen::MatrixXi Mesh(Ny,Nx); 
    discretize_domain(Mesh);


    Eigen::MatrixXd X(Nx,Ny); 
    Eigen::MatrixXd Y(Nx,Ny);
    meshgrid(hx,hy,X,Y); 

    //write meshgrid in a csv file
    std::cout <<"Trying to write to mesh.csv file the meshgrid values"<< std::endl;
    write_meshgrid_to_csv(X,Y);
    std::cout <<"--Mesh.csv filled with success!"<< std::endl;


    // RHS Vector
    Eigen::VectorXd RHS((Nx)*(Ny));
    // Boundary values
    std::cout<<"Dirichlet Boundary(case 0) or Neumann-Dirichlet Boundary(case 1)  or Neumann Boundary(case 2)? Press 0 or 1 or 2"<<std::endl;
    int case_h=0;
    std::cin>>case_h;
    Eigen::VectorXd U(Nx*Ny);
    Eigen::SparseMatrix<double> A(Nx*Ny, Nx*Ny);
    //Neumann
    switch (case_h){
        case 0 :
        {
        std::cout << "Enter the boundary values [u_x0, u_xL, u_y0, u_yL]: " << std::endl;
        double u_x0; // u(0,y)  
        double u_xL; // u(L,y)
        double u_y0; // u(x,0) 
        double u_yL; // u(x,L)
        std::cout<<"u_x0: ";std:: cin >> u_x0;
        std::cout<<"u_xL: "; std:: cin >> u_xL;
        std::cout<<"u_y0: " ;std:: cin >> u_y0;
        std::cout<<"u_yL: " ;std:: cin >> u_yL;
        dirichlett_boundary(u_x0,u_xL,u_y0,u_yL, RHS, Mesh);

    // System Matirx
        //  Eigen::SparseMatrix<double> A(Nx*Ny, Nx*Ny);
         matrix_assemply_dirichlet(A, Mesh);
    
    // Solve Linear System
        
        break;
        }
       
      
    
    //Neumann-Dirichlet
        case 1 :
        {
            std::cout << "For this case we take only one Neumann boundary which is for x=L " << std::endl;
        std::cout << "Enter the boundary values [u_x0, u_y0, u_yL]: " << std::endl;
        double u_x0; // u(0,y)  
            std::vector<double> g_xL; // u(L,y)
        double m;
        double u_y0; // u(x,0) 
        double u_yL; // u(x,L)
        std::cout<<"u_x0: ";std:: cin >> u_x0;
        std::cout<<"u_y0: " ;std::cin >> u_y0;
        std::cout<<"u_yL: " ;std:: cin >> u_yL;
        std::cout << "Enter Neumann boundary for x=L:" << std::endl;

        std::cout<<"du/dx(L,y): "<<std::endl;
        std::cin>>m;
        double u_xL;
        u_xL = hx*m;
        neumann_dirichlett_boundary(u_xL,u_x0,u_y0,u_yL, RHS, Mesh);

    // System Matirx
        // Eigen::SparseMatrix<double> A(Nx*Ny, Nx*Ny);
        matrix_assemply_neumann_dirichlet(A, Mesh);
        //std::cout << A << std::endl;
    
    // Solve Linear System
        
        // U=solve_system(A,RHS);
        break;
        }
        
    // Neumann-Boundary case
    //     case 2 :
    //     {
    //     std::cout << "For this case we consider all the boundaries as a Neumann boundary " << std::endl;
    //     std::cout << "Enter the boundary values [du/dy(x,0), du/dx(L,y),du/dx(y,0), du/dy(x,L)]: " << std::endl;
    //     double du_x0;   
    //     double du_y0; 
    //     double du_yL; 
    //     double du_xL;
    //     std::cout<<"0,y: ";std:: cin >> du_x0;
    //     std::cout<<"L,y: ";std:: cin >> du_xL;
    //     std::cout<<"x,0: " ;std:: cin >> du_y0;
    //     std::cout<<"x,L; " ;std:: cin >> du_yL;
    //    du_x0 = hx*du_x0;
    //    du_xL = hx*du_xL;
    //     du_y0 = hy*du_y0;
    //     du_yL = hy*du_yL;
    //     std::cout<<du_yL<<std::endl;
    //     neumann_dirichlett_boundary(du_xL,du_x0,du_y0,du_yL, RHS, Mesh);

    // // System Matirx
    //     // Eigen::SparseMatrix<double> A(Nx*Ny, Nx*Ny);
    //     matrix_assemply_neumann(A, Mesh);
    //     std::cout << A << std::endl;
    
    // Solve Linear System
        // U=solve_system(A,RHS);

        //}
    }

    std::cout << "Choose a numerical solver: (0/1/2)"<< std::endl;
    std::cout << "0) Eigen Solver"<<std::endl;
    std::cout << "1) Jacobi Method"<<std::endl;
    std::cout << "2) Gauss-Seidel Method"<<std::endl;

    int inp;
    std::cin >> inp;

    if (inp == 0)
    {
        std::cout << "You selected: Eigen Solver "<<std::endl;
        U = LU_method(A,RHS);
        Eigen::VectorXd b(U.size());
        b = A*U;

        double error;
        error = error_estimation_inf_norm(b,RHS);
        std::cout << "Eigen Error: "<< error << std::endl;
    }
    else if (inp == 1)
    {
        std::cout << "You selected: Jacobi Method "<<std::endl;
        U = jacobi_method(A,RHS);
        Eigen::VectorXd b(U.size());
        b = A*U;

        double error;
        error = error_estimation_inf_norm(b,RHS);
        std::cout << "Jacobi Error: "<< error << std::endl;
    }
    else if (inp == 2)
    {
        std::cout << "You selected: Gauss-Seidel Method "<<std::endl;
        U = gauss_seidel(A,RHS);
        Eigen::VectorXd b(U.size());
        b = A*U;

        double error;
        error = error_estimation_inf_norm(b,RHS);
        std::cout << "Gauss-Seidel Error: "<< error << std::endl;
    }
    
    
    // Save result to csv
    std::cout <<"Trying to write to solution.csv file the result vector values"<< std::endl;
    write_result_vector_to_csv(U);
    std::cout <<"--solution.csv filled with success!"<< std::endl; 
    std::cout <<"The solution is in the folder ./results/heat_equation_app_solutions"<< std::endl;
   
   /*
        The solver does not compute the error, the solver only computes a solution
        as per the bounday and mesh values entered by the user. 
        The error is computed in the test files with a known solution. 
    */
    return 0;
}