#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <iostream>
#include <vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<string>
#include <fstream>


class Linear_Solver
{   
    protected:
        Eigen::SparseMatrix<double> _A;
        Eigen::VectorXd _b;
        Eigen::VectorXd _result;

    public:

        Linear_Solver(Eigen::SparseMatrix<double> A,Eigen::VectorXd b);
        Eigen::VectorXd get_result();
        double error_estimation_inf_norm(Eigen::VectorXd &x_numerical,
                             Eigen::VectorXd &real_x) ;
        void write_result_vector_to_csv(std::string filename = "./results/heat_equation_app_solutions/numerical_solution.csv");

        virtual Eigen::VectorXd solve_system() = 0;
};



#endif

