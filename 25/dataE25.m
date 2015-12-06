a = 9;
b = 4;
c = 0;
if (spring_dof == top_dof +1)
    c = 0.01;
end


EA = 10;
l0 = sqrt(a*a + b*b + c*c);
k = (EA/l0)*(a/l0)^2;
k_spring = k*1.1;
alpha = a/l0;
ep = [EA, 1];
n1 = [0*b, 0*a, 0];
n2 = [b, a, c];
n3 = [2*b, a*0, 0];
%n_spring = n2;
coord = [n1; n2; n3];
%coord = [n1; n2; n3, n_spring];
coord = coord';
coord0 = coord;
%Enod = [1, 1, 2; 2, 2, 3; 3, 2, 4]; 
Enod = [1, 1, 2;
    2, 2, 3;];

Edof = [Enod(:,1), node_dof(Enod(:,2)), node_dof(Enod(:,3))];
nelm = 2;
nnod = length(coord);
ndof = nnod*3;
a=zeros(ndof,1);
top_dof = 5;
    
bc = [1 0;2 0;3 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0];

bc = [1 0;2 0;3 0;6 0;7 0;8 0;9 0];

if spring_dof == top_dof +1
    bc = [1 0; 2 0; 3 0; 7 0; 8 0; 9 0];
end



p1 = 2*alpha^3/(3*sqrt(3));
fmax = p1 * EA;
fmax = fmax*10;

nbr_steps = 100;

f0 = zeros(ndof, 1);
df = f0;
a = f0;
%a(top_dof+1)= c;

df(top_dof) = -fmax/nbr_steps;


[Ex0,Ey0,Ez0]=coordxtr(Edof,coord,node_dof((1:nnod)'),2);
%eldraw3(Ex0,Ey0,Ez0,[1 4 1]);
   
    
    