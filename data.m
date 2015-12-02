clear;
load('geom7e1.mat')
%edof : rad 1, x1, y1, x2, y2, x45, y45
%dof: xn vector place 2n-1
%dof: yn vector place 2n


E = ones(3,1)*210e9; %Är detta verkligen isotropiskt material?

eldraw2(ex,ey,[1 3 1]);