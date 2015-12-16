function J = jacobian(ec,xi,eta)

%J = jacobian(ec,xi,eta)
% ec = [x1 .. x4; y1 .. y4];

x=ec(1,:)';  y=ec(2,:)';

for i=1:4
    dNdxi(1,i) = dN_dxi(i,xi,eta);
    dNdeta(1,i) = dN_deta(i,xi,eta);
end

dx_dxi=dNdxi*x;
dx_deta=dNdeta*x;
dy_dxi=dNdxi*y;
dy_deta=dNdeta*y;

J = [dx_dxi, dx_deta;
    dy_dxi, dy_deta];
