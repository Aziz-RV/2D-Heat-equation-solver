#ifndef NEUMANNDIRICHLET_SOLVER_H
#define NEUMANNDIRICHLET_SOLVER_H

#include "pde_solver.h"
class Neumann_Dirichlet : public Pde_Solver
{
    private:
        double _u_x0; // u(0,y) 
        double _u_y0; // u(x,0)
        double _u_yL; // u(0,L)
        double _du_dxL; // du/dx at (x,y) = (L,y) // Neumann Boundary

    public:
        /* 
            Constructor for the Neumann_Dirichlet Problem, uses the constructor of the base class.
            Input: Domain values and meshsize (same as base class constructor)
        */
        Neumann_Dirichlet(const double x0, const double xL, const double y0,
                const double yL, const double hx, const double hy);


        /*
            Initializes the boundary values as per input of the user.
        */
        void set_boundary_values();


        /*
            Assembles the system matrix for the discretized Laplacian
            with Neumann-Dirichlett (mixed) boundary conditions.
        */
        void matrix_assembly() override;


        /*
            Assembles the right hand side vector (b) for the 
            Neumann-Dirichlett (mixed) boundary conditions.
        */
        void rhs_assembly() override;

};
     

#endif