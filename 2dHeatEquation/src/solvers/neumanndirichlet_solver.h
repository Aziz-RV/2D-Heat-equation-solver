#ifndef NEUMANNDIRICHLET_SOLVER_H
#define NEUMANNDIRICHLET_SOLVER_H

#include "pde_solver.h"
class Neumann_Dirichlet : public Pde_Solver
{
    private:
        double _u_x0;
        double _u_y0;
        double _u_yL;
        double _du_dxL;

    public:
        Neumann_Dirichlet(const double x0, const double xL, const double y0,
                const double yL, const double hx, const double hy);

        void set_boundary_values();
        void matrix_assembly() override;
        void rhs_assembly() override;

};
     

#endif