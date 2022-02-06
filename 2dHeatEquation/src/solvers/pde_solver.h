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
        std:: size_t _Nx; // Number of points in the mesh along x-axis
        std:: size_t _Ny; // Number of points in the mesh along y-axis
        std::size_t _N;
        const double _hx; 
        const double _hy;
        Eigen::MatrixXi _Mesh;
        Eigen::MatrixXd _X;
        Eigen::MatrixXd _Y;
        Eigen::SparseMatrix<double> _A;
        Eigen::VectorXd _RHS;

    public:

        /*
            Constructor for the Pde_Solver class. The constructor sets the rectangular domain, obatains the mesh, 
            prepares the meshgrid and initializes the system matrix _A and the _RHS vector.

            @param: x0 is the initial value of the domain along x-axis (initial x-coordinate).
            @param: y0 is the initial value of the domain along y-axis (initial y-coordinate).
            @param: xL is the final value of the domain along x-axis (final x-coordinate).
            @param: y0 is the initial value of the domain along y-axis (final y-coordinate).
            @param: hx is the meshsize along x-axis.
            @param: hy is the meshsize along y-axis..
        */
        Pde_Solver(const double x0, const double xL, const double y0,const double yL, const double hx, const double hy);
        
        
        /*
            Returns the sparse matrix _A.
        */
        Eigen::SparseMatrix<double>& get_A();


        /*
            Returns the vector _RHS.
        */
        Eigen::VectorXd& get_RHS();


        /*
            Discretizes the domain. (Assigns a number to each of the grid point 
            as per the lexicographical ordering).
        */
        void discretize_domain();

        
        /*
            Meshgrid is a set of two matrices _X and _Y which contain the x and y coordinates
            of the grid points respectively. Meshgrid is necessary when plotting the solution 
            of the pde using some scripting language like Matlab or Python.
        */
        void meshgrid();


        /*
            Generates a .csv file which contains the matrices _X and _Y (the meshgrid).
        */
        void write_meshgrid_to_csv(std::string filename = "./results/heat_equation_app_solutions/mesh.csv");


        /*
            Uses the linear solvers the compute solution to the linear system. Allows the user to choose
            from the LU solver (from Eigen), Jacobi solver and Gauss Seidel solver.
        */
        void render_solution();


        /*  
            Assembles the sparse system matrix. The actual assembly is dependent on the type of
            boundary conditions of the system, hence the actual implementation is done in derived classes. 
        */
        virtual void matrix_assembly() = 0;


        /*  
            Assembles the _RHS vector. The actual assembly is dependent on the type of
            boundary conditions of the system, hence the actual implementation is done in derived classes. 
        */
        virtual void rhs_assembly() = 0;
};



#endif // PDE_SOLVER_Η
