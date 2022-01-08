#ifndef LU_SOLVER_EIGEN_H
#define LU_SOLVER_EIGEN_H


#include "Linear_solver.h"


class LU_Solver_Eigen : public Linear_Solver
{
    public:
        LU_Solver_Eigen(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b);

        Eigen::VectorXd solve_system() override;

};

#endif