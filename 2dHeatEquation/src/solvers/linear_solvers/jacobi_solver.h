#ifndef JACOBI_SOLVER_H
#define JACOBI_SOLVER_H

#include "Linear_solver.h"

class Jacobi_Solver : public Linear_Solver
{ 
    public: 
        Jacobi_Solver(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b);
                    

        Eigen::VectorXd solve_system() override;


};

#endif