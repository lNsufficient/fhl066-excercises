function [ee, eff]  = plan3gs(ec,ed)
%Computes the Green-Lagrange strains and the deformation gradient
%   ee = [Exx; Eyy; 2Exy] enligt 5.24 (Matti)
%   eff = [dx/dx0; dx/dy0, dy/dx0; dy/dy0] 
%   ec - undeformed configuration [xrow; yrow]
%   ed - förskjutningar - a1 - x nod 1, a2 y nod 1, osv
evenIndex = (2:2:6);
oddIndex = (1:2:5);

x = ec(1,:) + ed(oddIndex);
y = ec(2,:) + ed(evenIndex);

A = abs((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)))/2;

dNdx = 1/(2*A)*[y(2)-y(3), y(3)-y(1), y(1)-y(2)];
dNdy = 1/(2*A)*[x(3)-x(2), x(1) - x(3), x(2) - x(1)];

B0 = zeros(3, 6);
B0([1;3], oddIndex) = [dNdx;dNdy];
B0([2;3], evenIndex) = [dNdy; dNdx];

H = zeros(4,6);
H([1:2]',[oddIndex]) = [dNdx;dNdy];
H([3:4]',[evenIndex]) = [dNdx;dNdy];

eff = H * ed;
A = [eff(1), 0, eff(3), 0;
    0, eff(2), 0, eff(4); 
    eff(2), eff(1), eff(4), eff(3)];
ee = (B0 +1/2*A*H)*ed;  
     
end

