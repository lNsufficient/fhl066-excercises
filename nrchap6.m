data
K0 = zeros(ndof, ndof);
f0 = zeros(ndof, 1);
an = zeros(ndof, 1);
fn = f0;
dbc = 10/n_end;
bc0 = ones(length(right_side), 1)*dbc;
bcn = bc; 

TOL = 1.2;
max_itr = 10;
for n=(1:n_end)
    bcn(right_side,2) = bcn(right_side,2) + bc0;
    a = an;
    %S = Sn;
    G = TOL + 1;
    j = 0;
    while(norm(G) > TOL && j < max_itr)
        j = j + 1;
        K = K0;
        fint = f0;
        for i = (1:nelm)
            dof = edof(i, 2:end);
            ed = a(dof);
            ec = [ex(i,:); ey(i,:)];
            [ee, ~] = plan3gs(ec, ed);
            es = D*ee;
            Ke = plan3ge(ec, t, D, ed, es);
            K(dof, dof) = K(dof, dof) + Ke;
            
            
            fe = plan3gf(ec, t, ed, es);
            fint(dof) = fint(dof) + fe;
        end
        G = fint - fn;
        da = solveq(K, -G, bc);
        a = a + da;
        %fint
        
    end
    an = a;
    %Sn = S;
end

plotpar= [1 1 2];
Ed = extract(edof, a);
eldisp2(ex, ey, Ed, plotpar);