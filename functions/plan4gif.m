function ef=plan4gif(ec,t,ed,es)

%ef=plan4gif(ec,t,ed,es)
% ---
% ef = [f_int1 .. f_int8];
% ec = [x1 .. x4; y1 .. y4]
% t tjocklek
% ed = [a1 .. a8];
% es = [S_xx; S_yy; S_xy]

gauss_points;

ef=zeros(1,8);


for i=1:4
    xi=xi_G(i);
    eta=eta_G(i);
    w=weight(i);
    [B0,A,H,~,J] = gradients4(ec,ed,xi,eta);
   
    B_u=A*H;
    
    B=B0+B_u;
    
    S=es{i};
    
    ef=ef+(B'*S*det(J)*t*w)';
    % ef är f_int'
end
