#include <iostream>
#include <vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<string>
#include <fstream>
#include "pde_solver.h"


        
Pde_Solver::Pde_Solver(const double x0, const double xL, const double y0,const double yL, const double hx, const double hy)
                : _x0(x0), _xL(xL), _y0(y0), _yL(yL), _hx(hx), _hy(hy)
        {
                
                    _Nx = (_xL - _x0) / _hx + 1;
                    _Ny = (_yL - _y0) / _hy + 1;
                    std::cout<<_Nx<<std::endl;
                    std::cout <<_Ny<<std::endl;
                    _N = _Nx * _Ny;
                    //initialize all Matrices
                    Eigen::MatrixXi Mesh(_Ny,_Nx); 
                    _Mesh=Mesh;
                    discretize_domain();
                    Eigen::MatrixXd X(_Nx,_Ny); 
                    Eigen::MatrixXd Y(_Nx,_Ny);
                    _X=X;
                    _Y=Y;
                    meshgrid();
                    Eigen::VectorXd RHS((_Nx)*(_Ny));
                    _RHS=RHS;
                    Eigen::SparseMatrix<double> A(_Nx*_Ny, _Nx*_Ny);
                    _A=A;
                    std::cout << "Setting up PDE Solver" << std::endl;
        }
        Eigen::SparseMatrix<double>& Pde_Solver::get_A(){
                return _A;}
        Eigen::VectorXd& Pde_Solver::get_RHS(){
                return _RHS;}
        void Pde_Solver::discretize_domain()
            {
                int m = 0;
                for (int i = 0; i < _Ny; i++)
                {
                    for (int j = 0; j < _Nx; j++)
                    {
                        _Mesh(i, j) = (i + 1) * m + (j + 1);
                    }
                    m++;
                }
            }
            

        void Pde_Solver::meshgrid()
            {
                for (int i = 0; i < _Ny; i++)
                {
                    for (int j = 0; j < _Nx; j++)
                    {
                        _X(i, j) = j * _hx;
                        _Y(i, j) = i * _hy;
                    }
                }
            }

        void Pde_Solver::write_meshgrid_to_csv(std::string filename)
            {   std::cout <<"Trying to write to mesh.csv file the meshgrid values"<< std::endl;
                std::ofstream csv_file(filename);

                for (int i = 0; i < _X.rows(); i++)
                {
                    for (int j = 0; j < _X.cols(); j++)
                    {
                        csv_file << _X(i, j) << ',' << _Y(i, j) << std::endl;
                    }
                }
                csv_file.close();
                std::cout <<"--Mesh.csv filled with success!"<< std::endl;
            }
        void Pde_Solver::render_solution(){
                int method=0;
                std::cout << "Choose a numerical solver: (0/1/2)"<< std::endl;
                std::cout << "0) Eigen Solver"<<std::endl;
                std::cout << "1) Jacobi Method"<<std::endl;
                std::cout << "2) Gauss-Seidel Method"<<std::endl;
                std::cin>>method;
                Eigen::VectorXd U;
                switch (method){
                    case 0:
                    {
                        LU_Solver_Eigen __LU_solver(_A,_RHS);
                        U = __LU_solver.solve_system();
                        Eigen::VectorXd b(U.size());
                        b = _A*U;
                        double error;
                        error = __LU_solver.error_estimation_inf_norm(b,_RHS);
                        std::cout << "Eigen Error: "<< error << std::endl;
                        __LU_solver.write_result_vector_to_csv();
                        break;
                    

                    }
                    case 1:
                    {
                        Jacobi_Solver __jacobi_solver(_A,_RHS);
                        U = __jacobi_solver.solve_system();
                        Eigen::VectorXd b(U.size());
                        b = _A*U;
                        double error;
                        error = __jacobi_solver.error_estimation_inf_norm(b,_RHS);
                        std::cout << "Jacobi Error: "<< error << std::endl;
                        __jacobi_solver.write_result_vector_to_csv();
                        break;
                    
                    }
                    case 2:
                    {
                        Gauss_Seidel_Solver __Gauss_Seidel_solver(_A,_RHS);
                        U = __Gauss_Seidel_solver.solve_system();
                        Eigen::VectorXd b(U.size());
                        b = _A*U;
                        double error;
                        error = __Gauss_Seidel_solver.error_estimation_inf_norm(b,_RHS);
                        std::cout << "Gauss Seidel Error: "<< error << std::endl; 
                        __Gauss_Seidel_solver.write_result_vector_to_csv();
                        break;
                        
                        
                    }
                }

                

            }




