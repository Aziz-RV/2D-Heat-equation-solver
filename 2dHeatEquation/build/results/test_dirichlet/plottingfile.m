clear;close
load('sol_1.csv')
load('mesh_1.csv')
mesh=mesh_1;
sol = sol_1;
X = mesh(:,1);
Y = mesh(:,2);
l = sqrt(length(sol));
sol = reshape(sol,[l,l]);
X = reshape(X,[l,l]);
Y = reshape(Y,[l,l]);
figure(1)
surf(X,Y,sol)
view(2)
shading interp 
colorbar
title('numeric')


f = sin(pi*X).*sinh(pi*Y);
figure(2)
surf(X,Y,f)
view(2)
shading interp
colorbar
title('real')

