data
K0 = zeros(ndof, ndof);
f0 = zeros(ndof, 1);
an = zeros(ndof, 1);
fn = f0;
bcn = bc;
bcn(right_side, 2) = 0;
maxItr = 10;
TOL = 1e-3;

for n=(1:n_end)

    a = an;
    %S = Sn;
    G = TOL + 1;
    bcn(right_side, 2) = dbc*n;
    size(bcn)
    a(right_side) = bcn(right_side, 2);
    size(a)
    j = 0;
    while(norm(G) > TOL && j < maxItr)
        j = j+1;
        K = K0;
        fint = f0;
        for i = (1:nelm)
            tmpedof = edof(i, 2:end);
            ed = a(tmpedof);
            ec = [ex(i,:); ey(i,:)];
            [ee, ~] = plan3gs(ec, ed);
            es = D*ee;
            Ke = plan3ge(ec, t, D, ed, es);
            K(tmpedof, tmpedof) = K(tmpedof, tmpedof) + Ke;
            
            fe = plan3gf(ec, t, ed, es);
            fint(tmpedof) = fint(tmpedof) + fe;
            max(abs(fe))
            pause;
            
        end
        G = -fint; % fint-0, fext = 0.

        G(right_side) = 0; %vi r�knr med att den blir bra
        %s�tt bc till att vara noll varje varv.
        da = solveq(K, G, bc);
        a = a + da;
        %fint
        
        G(bc(:,1)) = 0;
        bc(right_side, 2) = dbc*0;
    end
    if j == maxItr
        disp('reached maxitr')
    end

    an = a;
    %Sn = S;
end

plotpar= [1 1 2];
Ed = extract(tmpedof, a);
eldisp2(ex, ey, Ed, plotpar);