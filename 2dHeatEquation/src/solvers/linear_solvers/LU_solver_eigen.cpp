#include "LU_solver_eigen.h"


LU_Solver_Eigen::LU_Solver_Eigen(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b)
        :Linear_Solver(A,b)
    {
        std::cout<< "Preparing LU Solver" << std::endl;
    }


Eigen::VectorXd LU_Solver_Eigen::solve_system()
    {
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;
        start = std::chrono::system_clock::now();
        std::cout << "Solving system now" << std::endl;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>   solver;
        solver.analyzePattern(_A); 
        solver.factorize(_A); 
        // Eigen::VectorXd x = solver.solve(_b);
        // _result=x;
        // return x;
        _result = solver.solve(_b);
        end = std::chrono::system_clock::now();
        auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
         std::cout << "Elapsed time to solve the system: " << elapsed.count()/1000.0 << " seconds" << std::endl;
        return _result;
    }