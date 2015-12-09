function [es,ee] = bar3gs_updated(ec,ep,ed, old_eps)
%[es,ee] = bar3gs(ec,ep,ed)
% es normalkraft
% ee Greent�jning
% ---
% ec odeformerade koordinater [x1 x2; y1 y2 ; z1 z2]
% ep materialparametrar [E A0]
% ed f�rskjutningar [a1 a2 ... a6]
E=ep(1); A0 = ep(2);
L0 = norm(bar(ec));

%X = ec(:) + ed(:);
L = norm(bar_def(ec,ed));

de = (L^2-L0^2)/(2*L0^2);
ee = de + old_eps;

es=E*A0*ee;