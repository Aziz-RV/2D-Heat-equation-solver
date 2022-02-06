#include "solvers/pde_solver.h"
#include "solvers/neumanndirichlet_solver.h"
#include "solvers/dirichlet_solver.h"
#include "utilities.h"

/*
    The main file asks the user to enter the domain in which the pde will be solved.
    Then it prompts the user to enter the type of boundary conditions as well as the
    mesh-size that is needed.
*/

int main()
{
    conditions initial_conditions;
    initial_conditions =handle_user_inputs();

    // Boundary values
    std::cout<<"Dirichlet Boundary(case 0) or Neumann-Dirichlet Boundary(case 1), Press 0 or 1"<<std::endl;
    std::size_t case_h=0;
    std::cin>>case_h;
    switch(case_h){
        case 0 :
        {
        Dirichlet __PDE( initial_conditions.x0, initial_conditions.xL, initial_conditions.y0,initial_conditions.yL, initial_conditions.hx, initial_conditions.hy);
        __PDE.write_meshgrid_to_csv();
        __PDE.set_boundary_values();
        __PDE.rhs_assembly();
        __PDE.matrix_assembly();
        __PDE.render_solution();   
        break;
        
        }
        case 1:
        {
        Neumann_Dirichlet __PDE( initial_conditions.x0, initial_conditions.xL, initial_conditions.y0,initial_conditions.yL, initial_conditions.hx, initial_conditions.hy);
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
