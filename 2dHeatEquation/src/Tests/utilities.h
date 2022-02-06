#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<string>
#include <fstream>


void write_meshgrid_to_csv(Eigen::MatrixXd& X, Eigen::MatrixXd& Y,std::string filename = "./results/heat_equation_app_solutions/mesh.csv");
void write_result_vector_to_csv(Eigen::VectorXd& result, std::string filename = "./results/heat_equation_app_solutions/numerical_solution.csv");
void read_real_solution_from_csv(Eigen::VectorXd& _reference, std::string filename);
void input_user_mesh();
void discretize_domain(Eigen::MatrixXi& Mesh);
void meshgrid(double hx, double hy, Eigen::MatrixXd& X, Eigen::MatrixXd& Y);
void dirichlett_boundary(double u_x0, double u_xL, double u_y0, double u_yL ,Eigen::VectorXd& RHS, Eigen::MatrixXi& Mesh);
void neumann_dirichlett_boundary(double du_xL,double u_x0, double u_y0, double u_yL ,Eigen::VectorXd& RHS, Eigen::MatrixXi& Mesh);

void matrix_assemply_dirichlet( Eigen::SparseMatrix<double> &A, Eigen::MatrixXi& Mesh);
void matrix_assemply_neumann_dirichlet( Eigen::SparseMatrix<double> &A, Eigen::MatrixXi& Mesh);

double error_estimation(Eigen::VectorXd &x_numerical,Eigen::VectorXd &real_x );
double error_estimation_inf_norm(Eigen::VectorXd &x_numerical,Eigen::VectorXd &real_x );
Eigen::VectorXd LU_method(Eigen::SparseMatrix<double> &A,Eigen::VectorXd &vec);
Eigen::VectorXd gauss_seidel(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);
Eigen::VectorXd jacobi_method(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);
