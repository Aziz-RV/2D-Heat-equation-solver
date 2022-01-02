#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>


double error_estimation_inf_norm(Eigen::VectorXd &x_numerical,Eigen::VectorXd &real_x )
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

Eigen::VectorXd jacobi_method(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    int l = b.size();
    Eigen::VectorXd x(l);
    x = Eigen::VectorXd::Zero(l);

    Eigen::VectorXd y(l);
    y = Eigen::VectorXd::Zero(l);

    double tol = 0.00001;

    double ax = 0;
    Eigen::VectorXd bb(l);
    bb = A*x;
    int k = 0;

    while(error_estimation_inf_norm(bb,b)> tol && k < 1000)
    {
        for (int i = 0; i<l ;i++)
        {
            ax = 0;
            for (int j = 0; j < l ; j++)
            {   
                if(i != j)
                {
                    ax = ax + A.coeffRef(i,j)*x(j);
                }
            }
            y(i) = 1/A.coeffRef(i,i)*(b(i) - ax);
        }
        k = k+1;
        x = y;
        bb = A*x;
    }
    std::cout << "Total Iterates: " << k << std::endl;
    std::cout << "Result is: " << x << std::endl;
    return x;
}


Eigen::VectorXd gauss_seidel(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    int l = b.size();
    Eigen::VectorXd x(l);
    x = Eigen::VectorXd::Zero(l);

    Eigen::VectorXd y(l);
    y = Eigen::VectorXd::Zero(l);

    double tol = 0.00001;

    double ax = 0;
    double ay = 0;

    Eigen::VectorXd bb(l);
    bb = A*x;
    int k = 0;

    while(error_estimation_inf_norm(bb,b)> tol && k < 1000)
    {
        for (int i = 0; i<l ;i++)
        {
            ax = 0;
            ay = 0;
            for (int j = 0; j < i ; j++)
            {   
                ax = ax + A.coeffRef(i,j)*y(j);
            }
            for (int j = i+1; j<l ;j++)
            {
                ay = ay + A.coeffRef(i,j)*x(j);
            }
            y(i) = 1/A.coeffRef(i,i)*(b(i) - ax - ay);
        }
        k = k+1;
        x = y;
        bb = A*x;
    }
    std::cout << "Total Iterates: " << k << std::endl;
    std::cout << "Result is: " << x << std::endl;
    return x;
}


int main()
{
    Eigen::VectorXd vec(10);
    for(int i = 0; i < 10 ;i++)
    {
        vec(i) = 0.5*i;
    }
    std::cout << vec <<std::endl;

    Eigen::SparseMatrix<double> A(10,10);

    int l = vec.size();
    for(int i = 0; i < vec.size() ;i++)
    {
         A.coeffRef(i,i) = -2;
         if (i > 0 && i < l-1 )
         {
             A.coeffRef(i,i+1) = 1;
             A.coeffRef(i,i-1) = 1;
         }
    }
    A.coeffRef(0,1) = 1;
    A.coeffRef(l-1,l-2) = 1;

    std::cout << "Choose a numerical solver: (0/1/2)"<< std::endl;
    std::cout << "0) Eigen Solver"<<std::endl;
    std::cout << "1) Jacobi Method"<<std::endl;
    std::cout << "2) Gauss-Seidel Method"<<std::endl;

    int inp;
    std::cin >> inp;
    if (inp == 1)
    {
        Eigen::VectorXd x(l);
        x = jacobi_method(A,vec);
        Eigen::VectorXd b(l);
        b = A*x;

        double error;
        error = error_estimation_inf_norm(b,vec);
        std::cout << "Jacobi Error: "<< error << std::endl;
    }
    else if(inp == 2)
    {
        Eigen::VectorXd x(l);
        x = gauss_seidel(A,vec);
        Eigen::VectorXd b(l);
        b = A*x;

        double error;
        error = error_estimation_inf_norm(b,vec);
        std::cout << "GS Error: "<< error << std::endl;
    }
    return 0;
}

