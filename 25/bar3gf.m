function ef=bar3gf(ec,ed,es)
% ef ef=bar3gf(ec,ed,es) interna krafter i st�ngelementet 
% ---
% ec odeformerade koordinater [x1 x2; y1 y2 ; z1 z2]
% ed f�rskjutningar [a1 a2 ... a6]
% es normalkraft
L0=norm(bar(ec));

f_B = es/L0*bar_def(ec,ed);
f_A = -f_B;

ef = vertcat(f_A,f_B); %returnerar krafterna, sorterade i nod-ordning