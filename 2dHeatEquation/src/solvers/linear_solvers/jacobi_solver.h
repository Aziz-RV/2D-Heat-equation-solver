#ifndef JACOBI_SOLVER_H
#define JACOBI_SOLVER_H

#include "Linear_solver.h"

class Jacobi_Solver : public Linear_Solver
{ 
    public:

        /*
            Constructor for the Jacobi-Solver.
            Initializes the sparse matrix A and vector b in the system Ax = b
            using the constructor of base-class;
        */ 
        Jacobi_Solver(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b);
                    
        /*
            Uses the Jacobi Method to solve the sparse linear system.
            The solver has a tolerance of 0.000001 and the maximum numbers of interations is
            set by the user. The solver stops when either the error value is below the tolerance
            or when the iterations become equal to the maximum number of iterations. 
        */
        Eigen::VectorXd solve_system() override;


};

#endif