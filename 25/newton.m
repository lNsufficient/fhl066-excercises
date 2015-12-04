dataE25

f_n = f0;
f_l = f0;
K = zeros(ndof, ndof);

TOL = 1e-5;


plot_f = zeros(nbr_steps,1);
plot_u = zeros(nbr_steps,2);

k_spring = 0.001;
spring_dof = top_dof;
spring_vector = zeros(ndof,1);
spring_vector(spring_dof) = 1;

for n = [1: nbr_steps]
    f_l = f_l + df;
    f_n = f_l;
    %u_n = u_n;
    fint = df*0;
    r = TOL*norm(f_n) + 1;
    while(norm(r) > TOL*norm(f_n))
        
        K = 0*K;
        fint = 0*f_n;
        for i = [1:nelm]
            edof = Edof(i, :);
            edof = edof(2:end);
            ec = coord(edof);
            ec = [ec(1:3)', ec(4:6)'];
            ed = a(edof);
            [es, ~] = bar3gs(ec, ep, ed);
            Ke = bar3ge(ec,ep,ed,es);
            fe = bar3gf(ec, ed, es);
            K(edof, edof) = K(edof, edof) + Ke;
            fint(edof) = fint(edof) + fe;
        end
        K(top_dof, top_dof) = K(top_dof, top_dof) + [1 -1; -1 1]*k_spring;
        r = f_n - fint;
        r(bc(:,1)) = 0;
        da = solveq(K, r, bc);
        a = a + da;
    end
    plot_f(n) = f_eff(top_dof);
    plot_u(n,:) = a([top_dof, top_dof+1]);
end
coord0 = coord0';
coord = coord0 * 0;
for j=1:3
            coord(:,j) = coord0(:,j)+a(j:3:(end+j-3));
end

[Ex,Ey,Ez]=coordxtr(Edof,coord,node_dof((1:nnod)'),2);
%eldraw3(Ex,Ey,Ez,[1 4 1]);

plot3(plot_u(:,1), plot_u(:,2), (plot_f))
%plot3(plot_u(:,1), plot_u(:,1), (plot_f))
xlabel('förskjutning (Y) / meter')
ylabel('förskjutning (Z) / meter')
zlabel('kraft / Newton')
