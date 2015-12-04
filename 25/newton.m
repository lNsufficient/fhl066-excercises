dataE25

f_n = f0;
K = zeros(ndof, ndof);

TOL = 1e-5;

plot_f = zeros(nbr_steps,1);
plot_u = zeros(nbr_steps,1);

for n = [1: nbr_steps]
    f_n = f_n + df;
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
        r = f_n - fint;
        r(bc(:,1)) = 0;
        da = solveq(K, r, bc);
        a = a + da;
    end
    plot_f(n) = f_n(top_dof);
    plot_u(n) = a(top_dof);
end
coord0 = coord0';
coord = coord0 * 0;
for k=1:3
            coord(:,k) = coord0(:,k)+a(k:3:(end+k-3));
end

[Ex,Ey,Ez]=coordxtr(Edof,coord,node_dof((1:nnod)'),2);
%eldraw3(Ex,Ey,Ez,[1 4 1]);

plot((plot_u), (plot_f))
xlabel('förskjutning / meter')
ylabel('kraft / Newton')
