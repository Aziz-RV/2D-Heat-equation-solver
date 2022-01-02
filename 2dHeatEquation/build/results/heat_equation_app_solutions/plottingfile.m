clc;clear;close all;
load('numerical_solution.csv')
load('mesh.csv')
sol=numerical_solution;
X = mesh(:,1);
Y = mesh(:,2);
l = sqrt(length(sol));
sol = reshape(sol,[l,l]);
X = reshape(X,[l,l]);
Y = reshape(Y,[l,l]);
surf(X,Y,sol)
view(2)
shading interp
colorbar
