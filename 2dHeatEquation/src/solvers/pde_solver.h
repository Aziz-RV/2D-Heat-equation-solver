#ifndef PDE_SOLVER_H
#define PDE_SOLVER_H
#include <iostream>
#include <vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<string>
#include <fstream>
#include "linear_solvers/Linear_solver.h"
#include "linear_solvers/LU_solver_eigen.h"
#include "linear_solvers/jacobi_solver.h"
#include "linear_solvers/gauss_seidel_solver.h"


class Pde_Solver
{
    private:
        const double _x0;
        const double _xL;
        const double _y0;
        const double _yL;
        

    protected:
        int _Nx;
        int _Ny;
        int _N;
        const double _hx; // Needed in Neumann-Dirichlett
        const double _hy;
        Eigen::MatrixXi _Mesh;
        Eigen::MatrixXd _X;
        Eigen::MatrixXd _Y;
        Eigen::SparseMatrix<double> _A;
        Eigen::VectorXd _RHS;

    public:
        Pde_Solver(const double x0, const double xL, const double y0,const double yL, const double hx, const double hy);
        Eigen::SparseMatrix<double>& get_A();
        Eigen::VectorXd& get_RHS();
        void discretize_domain();
        void meshgrid();
        void write_meshgrid_to_csv(std::string filename = "./results/heat_equation_app_solutions/mesh.csv");
        void render_solution();
        virtual void matrix_assembly() = 0;
        virtual void rhs_assembly() = 0;
};



#endif // PDE_SOLVER_Η
