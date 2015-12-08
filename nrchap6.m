data
K0 = zeros(ndof, ndof);
f0 = zeros(ndof, 1);
for n=(1:n_end)
    fn = fn + df;
    a = an;
    S = Sn;
    G = TOL + 1;
    while(norm(G) < TOL)
        K = K0;
        fint = f0;
        for i = (1:nelm)
            edof = edof(i, 2:end);
            ed = a(edof);
            ec = [ex(k,:); ey(k,:)];
            [ee, ~] = plan3gs(ec, ed);
            es = D*ee;
            Ke = plan3ge(ec, t, D, ed, es);
            K(edof, edof) = K(edof, edof) + Ke;
            
            
            fe = plan3gf(ec, t, ed, es);
            fint(edof) = fint(edof) + fe;
        end
        da = K\(-G);
        a = a + da;
        %fint
        G = fint - fn;
    end
    an = a;
    Sn = S;
end