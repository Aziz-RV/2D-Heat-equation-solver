#ifndef DIRICHLET_SOLVER_H
#define DIRICHLET_SOLVER_H
#include "pde_solver.h"
class Dirichlet : public Pde_Solver
{
    private:
        double _u_x0;
        double _u_y0;
        double _u_yL;
        double _u_xL;

    public:
        Dirichlet(const double x0, const double xL, const double y0,
                const double yL, const double hx, const double hy);
        void set_boundary_values();
        void matrix_assembly() override;
        void rhs_assembly() override;
       

};

#endif