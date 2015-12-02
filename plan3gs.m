function [ee, eff]  = plan3gs(ec,ed)
%Computes the Green-Lagrange strains and the deformation gradient
%   ee = [Exx; Eyy; 2Exy] enligt 5.24 (Matti)
%   eff = [dx/dx0; dx/dy0, dy/dx0; dy/dy0] 
%   ec - undeformed configuration [xrow; yrow]
%   ed - förskjutningar - a1 - x nod 1, a2 y nod 1, osv
[B0, Au, H] = getB(ec, ed);
eff = H * ed;
ee = (B0 + Au*H)*ed;  
     
end

