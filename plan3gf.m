function ef = plan3gf(ec,t,ed,es)
%PLAN3GF ef = plan3gf(ec,t,ed,es) Computes internal element forces vector
%   ec - undeformed configuration [xrow; yrow]
%   ed - förskjutningar - a1 - x nod 1, a2 y nod 1, osv
%   es - S vector S11, S22, S33
%   t - thickness
%   ef - [f1, f2, ...,f6]
[x, y] = findxy(ec, ed)
vInt = calcArea(x, y)*t; 
[B0, Au, H]= getB(ec, ed);
B = B0 + 1/2*Au*H;
ef = B'*es*vInt; % B, es constant for the whole integral




end

