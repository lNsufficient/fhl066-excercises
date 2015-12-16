function Ke = plan3ge(ec,t,D,ed,es)
%Ke = plan3ge(ec,t,D,ed,es) provides element stiff matrix Ke for triag 3 nod large def element
%plane strain or plane stress
%   ec - undeformed configuration [xrow; yrow]
%   ed - förskjutningar - a1 - x nod 1, a2 y nod 1, osv
%   es - S vector S11, S22, S12
%   t - thickness
%   D - material properties
[B0, Au, H] = getB(ec, ed);
volIntegral0 = volumeInte(ec, zeros(numel(ec), 1), t);
S = [es(1), es(3); es(3), es(2)];
R = [S, zeros(size(S)); zeros(size(S)), S];
B = B0 + Au*H;
Ke = volIntegral0*(B'*D*B + H'*R*H);




end

