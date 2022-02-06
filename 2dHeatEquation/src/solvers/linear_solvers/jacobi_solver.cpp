#include "jacobi_solver.h"

Jacobi_Solver::Jacobi_Solver(Eigen::SparseMatrix<double> A,
                Eigen::VectorXd b)
    :Linear_Solver(A,b)
    {
        std::cout << "Preparing Jacobi Solver" << std::endl;
    }


Eigen::VectorXd Jacobi_Solver::solve_system()
{
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
    std::cout<< "max iterations :";
    std::size_t max_iter=0;
    std::cin >>max_iter;
    std::size_t l = _b.size();
    // Eigen::VectorXd x(l);
    // x = Eigen::VectorXd::Zero(l);
    _result = Eigen::VectorXd::Zero(l);

    Eigen::VectorXd y(l);
    y = Eigen::VectorXd::Zero(l);

    double tol = 0.000001;

    double ax = 0;
    Eigen::VectorXd bb(l);
    bb = _A*_result;
    int k = 0;

    while(error_estimation_inf_norm(bb,_b)> tol && k < max_iter)
    {
        for (int i = 0; i<l ;i++)
        {
            ax = 0;
            for (int j = 0; j < l ; j++)
            {   
                if(i != j)
                {
                    ax = ax + _A.coeffRef(i,j)*_result(j);
                }
            }
            y(i) = 1/_A.coeffRef(i,i)*(_b(i) - ax);
        }
        k = k+1;
        _result = y;
        bb = _A*_result;
    }
    end = std::chrono::system_clock::now();
     auto elapsed =
    std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Elapsed time to solve the system: " << elapsed.count()/1000.0 << " seconds" << std::endl;
    std::cout << "Total Iterates: " << k << std::endl;
    return _result;
    
}