#include "dirichlet_solver.h"
Dirichlet::Dirichlet(const double x0, const double xL, const double y0,
                const double yL, const double hx, const double hy):
        Pde_Solver(x0,  xL, y0,yL,  hx, hy)
{
    std::cout << "You have selected Dirichlet Boundaries"<< std::endl;
}

void Dirichlet::set_boundary_values()
{   
    std::cout << "Setting up Dirichlet Boundary values" << std::endl;
    std::cout << "Enter the boundary values [u(0,y), u(x,0), u(L,y), u(x,L)]: " << std::endl; 
    std::cout<<"u(0,y): ";std:: cin >> _u_x0;
    while(std::cin.fail()) {
           std::cout << "This is not a double, please enter a double for boundary values "<< std::endl;
           std::cin.clear();
           std::cin.ignore(256,'\n');
           std::cin >> _u_x0;
        }
    std::cout<<"u(x,0): " ;std::cin >> _u_y0;
     while(std::cin.fail()) {
           std::cout << "This is not a double, please enter a double for boundary values "<< std::endl;
           std::cin.clear();
           std::cin.ignore(256,'\n');
           std::cin >> _u_y0;
        }

    std::cout<<"u(L,y): " ;std:: cin >> _u_xL;
     while(std::cin.fail()) {
           std::cout << "This is not a double, please enter a double for boundary values "<< std::endl;
           std::cin.clear();
           std::cin.ignore(256,'\n');
           std::cin >> _u_xL;
        }
    std::cout<<"u(x,L): " ;std:: cin >> _u_yL;
     while(std::cin.fail()) {
           std::cout << "This is not a double, please enter a double for boundary values "<< std::endl;
           std::cin.clear();
           std::cin.ignore(256,'\n');
           std::cin >> _u_yL;;
        }
}


void Dirichlet::matrix_assembly() 
{
    std:: cout << "Started Matrix Assembly" << std::endl;

    int Lx = _Nx-1;
    _A.reserve(Eigen::VectorXi::Constant(_A.cols(),6));
    for (int j = 0; j< _Nx; j++) 
    {
        _A.coeffRef(j*_Ny,j*_Ny) = 1;
        _A.coeffRef(Lx + j*_Ny, Lx + j*_Ny) = 1;
    }

    int Ly = _Ny-1;
    for  (int i = 0; i < _Ny-1; i++)
    {
        _A.coeffRef(i,i) = 1;
        _A.coeffRef(Ly*_Nx + i,Ly*_Nx + i) = 1;
    }

    for (int i = 0; i < _N; i++)
    {
        if(_A.coeffRef(i,i) != 1)
        {
            _A.coeffRef(i,i) = -4;
            _A.coeffRef(i,i+1) = 1;
            _A.coeffRef(i,i-1) = 1;
            _A.coeffRef(i, i+_Ny) = 1;
            _A.coeffRef(i, i - _Ny) = 1;
        }
    }
     
}


void Dirichlet::rhs_assembly() 
{
    std:: cout << "Working on RHS vector" << std::endl;

    int x = 0; 
    int Lx = _Nx-1; 
    for (int j = 1; j < _Nx-1; j++) 
    {   
        _RHS(x + j*_Ny) = _u_x0;      
        _RHS(Lx + j*_Ny) = _u_xL;      
    } 
    
    int y = 0;
    int Ly = _Ny-1;
    for (int i = 0; i < _Ny; i++)
    {
        _RHS(y*_Nx + i) = _u_y0;
        _RHS(Ly*_Nx + i) = _u_yL;
    }
}
