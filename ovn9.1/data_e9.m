
load('task1_mesh2')
%Did not load data3....

%edof : rad 1, x1, y1, x2, y2, x45, y45
%dof: xn vector place 2n-1
%dof: yn vector place 2n

mm = 1e-3;
E = 210e3;
ex = ex;
ey = ey;

rho=7800e-9;

t = 1;
v = 0.3;

mp = [E, v];


ndof = nrdof;
nelm = nrelem;


%eldraw2(ex,ey,[1 3 1]);