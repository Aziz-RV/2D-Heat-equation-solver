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
        /*
            Constructor for Linear Solver. 
            Initializes the sparse matrix A and vector b in the system Ax = b.
        */
        Linear_Solver(Eigen::SparseMatrix<double> A,Eigen::VectorXd b);


        /*
            Returns the result vector.
        */
        Eigen::VectorXd get_result();


        /*
            Computes the error between two vectors with respect to the inf-norm.
            This function is used by Jacobi and Gauss_Seidel Solvers to compute error
            which is the compared to tolerance value.

            @ return: inf-norm difference of x and y.
        */
        double error_estimation_inf_norm(Eigen::VectorXd &x,
                             Eigen::VectorXd &y) ;


        /*
            Writes the result vector to .csv file
        */                     
        void write_result_vector_to_csv(std::string filename = "./results/heat_equation_app_solutions/numerical_solution.csv");

        /*
            Solves the sparse linear system to obtain result vector.
            This is a virtual function where the actual solver method is implemented in derived classes.

            @return: the solution (x) of linear system Ax = b.
        */
        virtual Eigen::VectorXd solve_system() = 0;
};



#endif

