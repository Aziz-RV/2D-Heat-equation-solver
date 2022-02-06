#include <iostream>
#include <vector>
#include<string>
#include <fstream>

//this structure the initial conditions given by the user for the Heat equation
struct conditions {
    double x0,xL,y0,yL,hx,hy;

};
/*
checks if the user typed a double
*/
void check_type_input(std::string msg, double _condition){
 while(std::cin.fail()) {
        std::cout << msg << std::endl;
        std::cin.clear();
        std::cin.ignore(256,'\n');
        std::cin >> _condition;
    }
}
/*
   This function checks the validity of user's input(the type of the input and also its coherence 
   with the conditions of the problem).
*/
conditions handle_user_inputs(){
    conditions initial_conditions;
    // x0 & xL
    std::cout << "Enter the domain:" << std::endl;
    std::cout << "Enter the initial and final x-value:" << std::endl;
    std::cout << "Enter the initial of x:" ;
    std::cin >> initial_conditions.x0;
    // checks if x0 typed is a double 
    check_type_input("This is not a double, please enter a double for intial x-value",initial_conditions.x0);
     
    std::cout << "and the final value of x:" ;
    // checks if xL typed is a double 
    std::cin >> initial_conditions.xL;
    check_type_input("This is not a double, please enter a double for final x-value",initial_conditions.xL);
    while (initial_conditions.x0 == initial_conditions.xL)
    {
        std::cout<< "Initial x-value must not equal final x-value."<<std::endl;
        std::cout << "Re-enter initial and final x-value:"<<std::endl;
        std::cout << "Enter the initial of x:" ;
        std::cin >> initial_conditions.x0;
        check_type_input("This is not a double, please enter a double for intial x-value",initial_conditions.x0);
        std::cout << "and the final value of x:" ;
        std::cin >> initial_conditions.xL; 
        check_type_input("This is not a double, please enter a double for final x-value",initial_conditions.xL);
      
    }

    //y0 && yL
    std::cout << "Enter the initial of y:" ;
    std::cin >> initial_conditions.y0;
    check_type_input("This is not a double, please enter a double for intial y-value",initial_conditions.y0);
    std::cout << "and the final value of y:" ;
    std::cin >> initial_conditions.yL;
    check_type_input("This is not a double, please enter a double for final y-value",initial_conditions.yL);
     while (initial_conditions.y0 == initial_conditions.yL)
    {
        std::cout<< "Initial x-value must not equal final y-value."<<std::endl;
        std::cout << "Re-enter initial and final y-value:"<<std::endl;
        std::cout << "Enter the initial of y:" ;
        std::cin >> initial_conditions.y0;
        check_type_input("This is not a double, please enter a double for intial y-value",initial_conditions.y0);
        std::cout << "and the final value of y:" ;
        std::cin >> initial_conditions.yL; 
        check_type_input("This is not a double, please enter a double for final y-value",initial_conditions.yL);
    }

    std::cout << "Enter the mesh size in x and y:" << std::endl;
    std::cout<<"hx= ";std::cin >> initial_conditions.hx;
    check_type_input("This is not a double, please enter a double for the x-Mesh-size",initial_conditions.hx);
   
    std::cout<<"hy= ";std::cin >> initial_conditions.hy;
    check_type_input("This is not a double, please enter a double for the y-Mesh-size",initial_conditions.hy);
    while(abs(initial_conditions.hx)>abs(initial_conditions.xL-initial_conditions.x0) || abs(initial_conditions.hy)> abs(initial_conditions.yL-initial_conditions.y0))
    {
        std::cout << "Error: Mesh-size in x is larger than domain size!"<<std::endl;
        std::cout << "Re-enter mesh-size"<< std::endl;
        std::cout<<"hx= ";std::cin >> initial_conditions.hx;
        check_type_input("This is not a double, please enter a double for the x-Mesh-size",initial_conditions.hx);
        std::cout<<"hy= ";std::cin >> initial_conditions.hy;
        check_type_input("This is not a double, please enter a double for the y-Mesh-size",initial_conditions.hy);
        
    }
    if (initial_conditions.x0>initial_conditions.xL) // if x0 > xL, mesh-size must be negative
    {
        initial_conditions.hx = -abs(initial_conditions.hx);
    }
    if (initial_conditions.y0>initial_conditions.yL)
    {
        initial_conditions.hy = -abs(initial_conditions.hy);
    }
    return initial_conditions;
}