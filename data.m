clear;
load('geom7e1.mat')
%edof : rad 1, x1, y1, x2, y2, x45, y45
%dof: xn vector place 2n-1
%dof: yn vector place 2n


E = 210e9; %Är detta verkligen isotropiskt material?
t = 1;
v = 0.3;
D = getD(E, v);
ndof = max(max(edof));
nelm = length(edof);
df = zeros(ndof, 1);


f_max = 20;
n_end = 200;
right_side = find(bc(:,2)> 0);
bc(right_side, 2) = 0; 


%eldraw2(ex,ey,[1 3 1]);