#ifndef GAUSS_SEIDEL_SOLVER_H
#define GAUSS_SEIDEL_SOLVER_H

#include "Linear_solver.h"

class Gauss_Seidel_Solver : public Linear_Solver
{
    public:
        Gauss_Seidel_Solver(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b);

        Eigen::VectorXd solve_system() override;

};

#endif