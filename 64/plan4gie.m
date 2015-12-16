function Ke=plan4gie(ec,t,D,ed,es)

%Ke=plan4gie(ec,t,D,ed,es)

a_x=ed(1,1:2:(end-1))';
a_y=ed(1,2:2:end)';

gauss_points;

K=zeros(8);

for i=1:4   %varje i motsvarar en gausspunkt
    xi=xi_G(i);
    eta=eta_G(i);
    w=weight(i);
    [B0,A,H,~,J] = gradients4(ec,ed,xi,eta);
    
    B_u = A*H;

    B = B0+B_u;
    
    es_point = es{i};
    D_point = D{i};
    
    s=es{i};
    if numel(s)==1
        keyboard
    end
    S = [s(1), s(3); s(3), s(2)];
    
    R=[S, eye(2); eye(2), S];
    
    K_incr = (B'*D_point*B+H'*R*H)*t*det(J)*w;
    K = K+K_incr;
end

Ke=K; 
