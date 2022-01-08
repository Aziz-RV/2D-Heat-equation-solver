#include "solvers/pde_solver.h"
#include "solvers/neumanndirichlet_solver.h"
#include "solvers/dirichlet_solver.h"

// we can add all the includes but they are already in utilities.h


int main()
{
    std::cout << "Enter the domain:" << std::endl;
    std::cout << "Enter the initial and final x-value:" << std::endl;
    const double x0;
    const double xL;
    std::cout << "Enter the initial of x:" ;
    std::cin >> x0;
    std::cout << "and the final value of x:" ;
    std::cin >> xL;
    const double y0;
    const double yL;
    std::cout << "Enter the initial of y:" ;
    std::cin >> y0;
    std::cout << "and the final value of y:" ;
    std::cin >> yL;
    const double hx; // mesh size
    const double hy; // mesh size
    std::cout << "Enter the mesh size in x and y:" << std::endl;
    std::cout<<"hx= ";std::cin >> hx;
    std::cout<<"hy= ";std::cin >> hy;
    int Nx = (xL-x0)/hx +1; // No. of columns (must be int)
    int Ny = (yL-y0)/hy +1; // No. of rows

    // Boundary values
    std::cout<<"Dirichlet Boundary(case 0) or Neumann-Dirichlet Boundary(case 1), Press 0 or 1"<<std::endl;
    int case_h=0;
    std::cin>>case_h;
    switch(case_h){
        case 0 :
        {
        Dirichlet __PDE( x0, xL, y0,yL, hx, hy);
        __PDE.write_meshgrid_to_csv();
        __PDE.set_boundary_values();
        __PDE.rhs_assembly();
        __PDE.matrix_assembly();
        __PDE.render_solution();   
        break;
        
        }
        case 1:
        {
        Neumann_Dirichlet __PDE( x0, xL, y0,yL, hx, hy);
        __PDE.write_meshgrid_to_csv();
        __PDE.set_boundary_values();
        __PDE.rhs_assembly();
        __PDE.matrix_assembly();
        __PDE.render_solution();
        break; 
        }

       
    
    }
    return 0;
}