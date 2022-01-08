#include "gauss_seidel_solver.h"

Gauss_Seidel_Solver::Gauss_Seidel_Solver(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b)
                :Linear_Solver(A,b)
    {
        std::cout << "Preparing Gauss-Seidel Solver"<<std::endl;
    }

Eigen::VectorXd Gauss_Seidel_Solver:: solve_system() 
{
    std::cout<< "max iterations :";
    int max_iter=0;
    std::cin >>max_iter;
    int l = _b.size();
    Eigen::VectorXd x(l);
    x = Eigen::VectorXd::Zero(l);

    Eigen::VectorXd y(l);
    y = Eigen::VectorXd::Zero(l);

    double tol = 0.000001;

    double ax = 0;
    double ay = 0;

    Eigen::VectorXd bb(l);
    bb = _A*x;
    int k = 0;

    while(error_estimation_inf_norm(bb,_b)> tol && k < max_iter)
    {
        for (int i = 0; i<l ;i++)
        {
            ax = 0;
            ay = 0;
            for (int j = 0; j < i ; j++)
            {   
                ax = ax + _A.coeffRef(i,j)*y(j);
            }
            for (int j = i+1; j<l ;j++)
            {
                ay = ay + _A.coeffRef(i,j)*x(j);
            }
            y(i) = 1/_A.coeffRef(i,i)*(_b(i) - ax - ay);
        }
        k = k+1;
        x = y;
        bb = _A*x;
    }
    std::cout << "Total Iterates: " << k << std::endl;
    // std::cout << "Result is: " << x << std::endl;
    _result=x;
    return x;
}
