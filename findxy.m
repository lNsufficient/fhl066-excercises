function [x, y] = findxy(ec, ed)
%FINDXY return x and y coordinates
%   ec - undeformed configuration [xrow; yrow]
%   ed - förskjutningar - a1 - x nod 1, a2 y nod 1, osv

x = ec(1,:) + ed(1:2:5)';
y = ec(2,:) + ed(2:2:6)';

end

