#include "Linear_solver.h"

Linear_Solver::Linear_Solver(Eigen::SparseMatrix<double> A,
                    Eigen::VectorXd b)
    {
        _A = A;
        _b = b;
            
    }

Eigen::VectorXd Linear_Solver::get_result()
    {
        return _result;
    }

double Linear_Solver::error_estimation_inf_norm(Eigen::VectorXd &x_numerical,
                                Eigen::VectorXd &real_x) 
    {
        double error = 0;
        double num;
        for (int i = 0; i<real_x.size(); i++)
        {
            num = abs(x_numerical(i)- real_x(i));
            if (num > error)
            {
                error = num;
            }
        } 
        return error;
    }


void Linear_Solver::write_result_vector_to_csv(std::string filename)

    {

        std::cout<<_result<<std::endl;
        std::cout <<"Trying to write to solution.csv file the result vector values"<< std::endl;
            
        std::ofstream csv_file(filename);
        for (int i = 0; i < _result.size() ; i++)
        {
            csv_file << _result(i) << std::endl;
        }
        csv_file.close();
        std::cout <<"--solution.csv filled with success!"<< std::endl; 
        std::cout <<"The solution is in the folder ./results/heat_equation_app_solutions"<< std::endl;
    };