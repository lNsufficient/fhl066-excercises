function [B0,A,H,displ_grad,J] = gradients4(ec,ed,xi,eta)
% Beräknar diverse derivataliknande storheter i punkten (xi,eta).
% [B0,A,H,displ_grad,J] = gradients4(ec,ed,xi,eta)
% displ_grad = [du_x_dx; du_x_dy; du_y_dx; du_y_dy]

a_x=ed(1,1:2:(end-1))';
a_y=ed(1,2:2:end)';

B0 = [];
H=[];
J = jacobian(ec,xi,eta);
J_inv = inv(J)';
for j=1:4   % varje j motsvarar en nod m. formfunkt.
    N_grad = J_inv*[dN_dxi(j,xi,eta);dN_deta(j,xi,eta)];
    dN_dx = N_grad(1);  dN_dy=N_grad(2);
    dN_dx_tot(1,j) = dN_dx;
    dN_dy_tot(1,j) = dN_dy;
    B_incr=[dN_dx,  0;
            0,  dN_dy;
            dN_dy, dN_dx];
    B0 = horzcat(B0,B_incr);
    H=horzcat(H,kron([1 0; 0 1], [dN_dx;dN_dy]));
    shapegrad_xi(1,j) = dN_dxi(j,xi,eta);
    shapegrad_eta(1,j) = dN_deta(j,xi,eta);
end

du_x_dx = dN_dx_tot*a_x;
du_y_dx = dN_dx_tot*a_y;
du_x_dy = dN_dy_tot*a_x;
du_y_dy = dN_dy_tot*a_y;

displ_grad = [du_x_dx; du_x_dy; du_y_dx; du_y_dy];

A = [du_x_dx,   0   , du_y_dx,  0   ;
    0   ,du_x_dy,   0   ,du_y_dy;
    du_x_dy,du_x_dx, du_y_dy, du_y_dx];
