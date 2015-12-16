function volumeIntegral = volumeInte(ec, ed, t)
%function volumeIntegral = volumeInte(ec, ed, t)
%   ec - undeformed configuration [xrow; yrow]
%   ed - förskjutningar - a1 - x nod 1, a2 y nod 1, osv
%   t - thickness
[x, y] = findxy_def(ec, ed);
volumeIntegral = calcArea(x, y)*t;

end

