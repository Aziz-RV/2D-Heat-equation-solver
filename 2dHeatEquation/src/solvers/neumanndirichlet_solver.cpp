#include "neumanndirichlet_solver.h"
Neumann_Dirichlet::Neumann_Dirichlet(const double x0, const double xL, const double y0,
                const double yL, const double hx, const double hy):
        Pde_Solver(x0,  xL, y0,yL,  hx, hy){};


void Neumann_Dirichlet::set_boundary_values()
        {
            std::cout << "Setting up Neumann-Dirichlet Boundary values" << std::endl;
            std::cout << "For this case we take only one Neumann boundary which is for x=L " << std::endl;
            std::cout << "Enter the boundary values [u_x0, u_y0, u_yL]: " << std::endl; 
             double m;
            std::cout<<"u_x0: ";std:: cin >> _u_x0;
            std::cout<<"u_y0: " ;std::cin >> _u_y0;
            std::cout<<"u_yL: " ;std:: cin >> _u_yL;
            std::cout << "Enter Neumann boundary for x=L:" << std::endl;
            std::cout<<"du/dx(L,y): "<<std::endl;
            std::cin>>m;
             double u_xL;
            _du_dxL = _hx*m;
        }

        void Neumann_Dirichlet::matrix_assembly() 
        {
            std:: cout << "Started Matrix Assembly" << std::endl;
            /*
                This function constructs system matrix for the discretized system.
                This system matrix takes depends on the boundary conditions and is the discretized form
                of the laplacian operator.
            */
            int Ny = _Mesh.rows();
            int Nx = _Mesh.cols();
            int N = (Nx)*(Ny);

            int Lx = Nx-1;
            _A.reserve(Eigen::VectorXi::Constant(_A.cols(),6));
            for (int j = 1; j< Nx-1; j++)
            {
                _A.coeffRef(j*Ny,j*Ny) = 1;
                _A.coeffRef(Lx + j*Ny, Lx + j*Ny) = 1;
                _A.coeffRef(Lx + j*Ny, (Lx + j*Ny)-1) = -1;
            }

            int Ly = Ny-1;
            for  (int i = 0; i < Ny; i++)
            {
                _A.coeffRef(i,i) = 1;
                _A.coeffRef(Ly*Nx + i,Ly*Nx + i) = 1;
            }

            for (int i = 0; i < N; i++)
            {
                if(_A.coeffRef(i,i) != 1)
                {
                    _A.coeffRef(i,i) = -4;
                    _A.coeffRef(i,i+1) = 1;
                    _A.coeffRef(i,i-1) = 1;
                    _A.coeffRef(i, i+ Ny) = 1;
                    _A.coeffRef(i, i - Ny) = 1;
                }
            }
        }

void Neumann_Dirichlet::rhs_assembly() 
        {
            std:: cout << "Working on RHS vector" << std::endl;
            int Nx = _Mesh.cols();
            int Ny = _Mesh.rows();

            int x = 0; // First Column
            int Lx = Nx-1; // Last column
            for (int j = 1; j < Nx-1; j++) // changed things here
            {   
                _RHS(x + j*Ny) = _u_x0;       // Makes First column 0 // u(0,y) =  0
                //std::cout<< x + j*Ny << ':' << u_x0 << std::endl;
                _RHS(Lx + j*Ny) = _du_dxL;      // Makes Last Column L // (du/dx)(L,y) = g(j)
                //std::cout<< Lx + j*Ny << ':' << u_xL << std::endl;
            } 
            
            int y = 0;
            int Ly = Ny-1;
            for (int i = 0; i < Ny; i++)
            {
                _RHS(y*Nx + i) = _u_y0;
                //std::cout<< y*Nx + i << ':' << u_y0 << std::endl;
               _RHS(Ly*Nx + i) = _u_yL;
                //std::cout<< Ly*Nx + i << ':' << u_yL << std::endl;
            }
        }

