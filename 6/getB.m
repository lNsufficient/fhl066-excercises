function [B0, Au, H] = getB(ec, ed)
%GETB Finds the B = getB(ec, ed)
%   ec - undeformed configuration [xrow; yrow]
%   ed - förskjutningar - a1 - x nod 1, a2 y nod 1, osv

evenIndex = (2:2:6);
oddIndex = (1:2:5);

%[x_def, y_def] = findxy_def(ec, ed);
[x0, y0] = findxy0(ec);
A = calcArea(x0, y0);

dNdx = 1/(2*A)*[y0(2)-y0(3), y0(3)-y0(1), y0(1)-y0(2)];
dNdy = 1/(2*A)*[x0(3)-x0(2), x0(1)-x0(3), x0(2)-x0(1)];

B0 = zeros(3, 6);
B0([1;3], oddIndex) = [dNdx;dNdy];
B0([2;3], evenIndex) = [dNdy; dNdx];

H = zeros(4,6);
H([1:2]',[oddIndex]) = [dNdx;dNdy];
H([3:4]',[evenIndex]) = [dNdx;dNdy];

eff = H * ed;
Au = [eff(1), 0, eff(3), 0;
    0, eff(2), 0, eff(4); 
    eff(2), eff(1), eff(4), eff(3)];


end

