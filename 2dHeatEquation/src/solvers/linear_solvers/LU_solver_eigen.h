#ifndef LU_SOLVER_EIGEN_H
#define LU_SOLVER_EIGEN_H


#include "Linear_solver.h"


class LU_Solver_Eigen : public Linear_Solver
{
    public:

        /*
            Constructor for the LU_Eigen-Solver.
            Initializes the sparse matrix A and vector b in the system Ax = b
            using the constructor of base-class;
        */ 
        LU_Solver_Eigen(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b);

        /*
            Uses the already implemented LU solver from the Eigen library
            to solve the linear system.
        */            
        Eigen::VectorXd solve_system() override;

};

#endif