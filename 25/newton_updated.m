top_dof = 5;
spring_dof = top_dof;

dataE25

f_n = f0;
f_l = f0;
K = zeros(ndof, ndof);

TOL = 1e-4;

%k_spring = 0;
plot_f = zeros(nbr_steps,1);
plot_u = zeros(nbr_steps,3);




if (spring_dof == 6)
    k_spring = 0.3;
end

k_spring = 0;
eps_n = 0;
eps_spring = 0;

for n = [1: nbr_steps]
    f_l = f_l + df;
    f_n = f_l;
    du_n = 0*a;

    r = TOL*norm(f_n) + 1;
    while(norm(r) > TOL*norm(f_l))
        
        K = 0*K;
        fint = 0*f_n;
        for i = [1:nelm]
            edof = Edof(i, :);
            edof = edof(2:end);
            ec = coord(edof);
            ec = [ec(1:3)', ec(4:6)'];
            ed = du_n(edof);
            [es, eps_n] = bar3gs_updated(ec, ep, ed, eps_n);
            Ke = bar3ge_updated(ec,ep,ed,es);
            fe = bar3gf(ec, ed, es);
            K(edof, edof) = K(edof, edof) + Ke;
            fint(edof) = fint(edof) + fe;
        end
        if spring_dof ~= -1
            K(spring_dof, spring_dof) = K(spring_dof, spring_dof) + k_spring;
            eps_spring = k_spring*du_n(spring_dof) + eps_spring;
            f_n(spring_dof) = f_l(spring_dof) - eps_spring;
        end
        r = f_n - fint
        da = solveq(K, r, bc)
        r(bc(:,1)) = 0;
        du_n = du_n + da;
        
        coord(:) = coord(:) + du_n(:);
        
    end
    a = a + du_n;
    plot_f(n) = f_n(top_dof);
    plot_u(n,:) = a([top_dof-1, top_dof, top_dof+1]);
end

[Ex,Ey,Ez]=coordxtr(Edof,coord,node_dof((1:nnod)'),2);
%eldraw3(Ex,Ey,Ez,[1 4 1]);
if spring_dof == top_dof+1
    plot3(plot_u(:,2), plot_u(:,3), (plot_f))
    xlabel('förskjutning (Y) / meter')
    ylabel('förskjutning (Z) / meter')
    zlabel('kraft / Newton')
elseif spring_dof == top_dof
    plot(plot_u(:,2), plot_f);
    xlabel('förskjutning (Y) / meter')
    ylabel('kraft / Newton')
    
elseif spring_dof == -1
    plot3(plot_u(:,1), plot_u(:,2), plot_u(:,3))
    xlabel('förskjutning (X) / meter')
    ylabel('förskjutning (Y) / meter')
    zlabel('förskjutning (Z) / meter')
end
%plot3(plot_u(:,1), plot_u(:,1), (plot_f))

