a = 9;
b = 4;
c = 0;

EA = 10;
l0 = sqrt(a*a + b*b + c*c);
alpha = a/l0;
ep = [EA, 1];
n1 = [0*b, 0*a, c];
n2 = [b, a, c];
n3 = [2*b, a*0, c];
coord = [n1; n2; n3];
coord = coord';
coord0 = coord;
Enod = [1, 1, 2;
    2, 2, 3];
Edof = [Enod(:,1), node_dof(Enod(:,2)), node_dof(Enod(:,3))];
nelm = 2;
nnod = length(coord);
ndof = nnod*3;
a=zeros(ndof,1);
top_dof = 5;
    
bc = [1 0;
    2 0;
    3 0;
    6 0; %Den kan inte påverkas i z-led i nuläget, får inte räkna med z.
    7 0;
    8 0;
    9 0];






p1 = 2*alpha^3/(3*sqrt(3));
fmax = p1 * EA;

nbr_steps = 50;

f0 = zeros(ndof, 1);
df = f0;
a = f0;

df(top_dof) = -fmax/nbr_steps;


[Ex0,Ey0,Ez0]=coordxtr(Edof,coord,node_dof((1:nnod)'),2);
%eldraw3(Ex0,Ey0,Ez0,[1 4 1]);
   
    
    