clear;close
load('sol_4.csv')
load('mesh_4.csv')
mesh=mesh_4;
sol = sol_4;
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


f = cos(pi*X).*sinh(pi*Y);
figure(2)
surf(X,Y,f)
view(2)
shading interp
colorbar
title('real')

