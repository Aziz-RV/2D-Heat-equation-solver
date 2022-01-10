#ifndef DIRICHLET_SOLVER_H
#define DIRICHLET_SOLVER_H
#include "pde_solver.h"
class Dirichlet : public Pde_Solver
{
    private:
        double _u_x0; // u(0,y) 
        double _u_y0; // u(x,0)
        double _u_yL; // u(0,L)
        double _u_xL; // u(L,y)

    public:
        /* 
            Constructor for the Dirichlet boundary conditions, uses the constructor of the base class.
            Input: Domain values and meshsize (same as base class constructor)
        */
        Dirichlet(const double x0, const double xL, const double y0,
                const double yL, const double hx, const double hy);


        /*
            Initializes the boundary values as per input of the user.
        */        
        void set_boundary_values();


        /*
            Assembles the system matrix for the discretized Laplacian
            with Dirichlett boundary conditions.
        */        
        void matrix_assembly() override;


        /*
            Assembles the right hand side vector (b) for the 
            Dirichlett boundary conditions.
        */        
        void rhs_assembly() override;
       

};

#endif