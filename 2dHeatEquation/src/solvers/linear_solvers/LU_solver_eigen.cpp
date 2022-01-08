#include "LU_solver_eigen.h"


LU_Solver_Eigen::LU_Solver_Eigen(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b)
        :Linear_Solver(A,b)
    {
        std::cout<< "Preparing LU Solver" << std::endl;
    }


Eigen::VectorXd LU_Solver_Eigen::solve_system()
    {
        std::cout << "Solving system now" << std::endl;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>   solver;
        solver.analyzePattern(_A); 
        solver.factorize(_A); 
        Eigen::VectorXd x = solver.solve(_b);
        _result=x;
        return x;
    }