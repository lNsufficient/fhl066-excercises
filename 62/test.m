load('control.mat')
ec = ec';
ed = ed';
[ee, eff] = plan3gs(ec, ed);
E_xx=ee(1); E_yy=ee(2); E_xy=ee(3)/2;

E=[E_xx, E_xy, 0;
   E_xy, E_yy, 0;
   0,0,0];

F=F_vect2tens(eff);

C=F'*F;


E
(E-(C-eye(3))/2)