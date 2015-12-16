function [ee,eff] = plan4gis(ec,ed)

% [ee,eff] = plan4gis(ec,ed)
% ee = {ee{1} .. ee{4}}
%   - ee{i} = [E_xx; E_yy; 2*E_xy] Greentöjningar per Gausspunkt
% eff = {eff{1} .. eff{4}}
%   - eff{i} = [dx_dx0; dx_dy0; dy_dx0; dy_dy0] deformationsgradient

gauss_points;

ee=cell(4,1);
eff=cell(4,1);

for i=1:4
    xi=xi_G(i);
    eta=eta_G(i);
    [B0,A,H,displ_grad,J] = gradients4(ec,ed,xi,eta);
    B_u=A*H;
    ee{i} = (B0 + B_u/2)*ed';
    eff{i} = displ_grad + [1;0;0;1];
end
