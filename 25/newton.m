top_dof = 5;
spring_dof = top_dof;

dataE25

f_n = f0;
f_l = f0;
K = zeros(ndof, ndof);

TOL = 1e-4;

%k_spring = 0;
k_spring = k_spring*4
plot_f = zeros(nbr_steps,1);
plot_u = zeros(nbr_steps,3);


%k_spring = 0;

use_hookes = 0;

coordprint = zeros(size(coord0))
if (spring_dof == 6)
    k_spring = 0.3;
end
    
for n = [1: nbr_steps]
    f_l = f_l + df;
    f_n = f_l;
    %u_n = u_n;

    r = TOL*norm(f_n) + 1;
    while(norm(r) > TOL*norm(f_l))
        
        K = 0*K;
        fint = 0*f_n;
        for i = [1:nelm]
            edof = Edof(i, :);
            edof = edof(2:end);
            ec = coord0(edof);
            ec = [ec(1:3)', ec(4:6)'];
            ed = a(edof);
            if use_hookes           
                [es, ~] = bar3gs(ec, ep, ed);
                Ke = bar3ge(ec,ep,ed,es);
                fe = bar3gf(ec, ed, es);
            else
                %[es, ee] = bar3gs_log1(ec, ep, ed)
                [es, ee] = bar3gs_log(ec, ep, ed)
                %Ke = bar3ge_log1(ec,ep,ed,es,ee)
                Ke = bar3ge_log(ec,ep,ed,ee)
                fe = bar3gf(ec, ed, es)
            end
            K(edof, edof) = K(edof, edof) + Ke;
            fint(edof) = fint(edof) + fe;
        end
        
        if spring_dof == 6
            K(spring_dof, spring_dof) = K(spring_dof, spring_dof) + k_spring/10;
            f_n(spring_dof) = f_l(spring_dof) - k_spring/10*a(spring_dof);
            
            K(top_dof, top_dof) = K(top_dof, top_dof) + k_spring;
            f_n(top_dof) = f_l(top_dof) - k_spring*a(top_dof);
        elseif spring_dof ~= -1
            K(spring_dof, spring_dof) = K(spring_dof, spring_dof) + k_spring;
            f_n(spring_dof) = f_l(spring_dof) - k_spring*a(spring_dof);
        end
        r = f_n - fint;
        da = solveq(K, r, bc);
        r(bc(:,1)) = 0;
        a = a + da;
        
        %for k=1:3
        %    coord(:,k) = coord0(:,k)+a(k:3:(end+k-3)); %Detta är NOG fel, men blir rätt för nod 5.
        %end
    end
    plot_f(n) = f_n(top_dof);
    plot_u(n,:) = a([top_dof-1, top_dof, top_dof+1]);
    
    coordprint(:) = coord0(:) + a(:);
    [Ex,Ey,Ez]=coordxtr(Edof,coordprint',node_dof((1:nnod)'),2);
    %eldraw3(Ex,Ey,Ez,[1 4 1]);
    %pause;
end

coord0 = coord0';
coord = coord0 * 0;
for j=1:3
            coord(:,j) = coord0(:,j)+a(j:3:(end+j-3));
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

