function f = spring_force(s0, x0, u, k)
%f = spring_force(s0, x0, u, k), calculates f_eff when a spring is attached from s0 to node at x
%   x0 nodes starting position
%   u  nodes displacement
%   s0 springs other coordinate
%   k spring constant

v0 = x0 - s0;

x = x0 + u;
v = x - s0;

l0 = sqrt(v0.*v0);
l = sqrt(x.*x);
dl = l - l0;
force = dl * k;
f = v*force/norm(force);


end

