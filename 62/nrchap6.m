data
K0 = zeros(ndof, ndof);
f0 = zeros(ndof, 1);
an = zeros(ndof, 1);

bcn = bc;
bcn(right_side, 2) = 0;
maxItr = 10;
TOL = 1e-6;

fplot = zeros(n_end, numel(right_side));
uplot = fplot*0;
gplot = [];
for n=(1:n_end)
    if mod(n, 2)
        disp(n)
    end
    disp('==========')
    a = an;
    %S = Sn;
    G = TOL + 1;
    bcn(right_side, 2) = dbc;
    %a(right_side_nodes) = bcn(right_side, 2);
    j = 0;
    while(norm(G) > TOL && j < maxItr || j < 2)
        j = j+1;
        K = K0*0;
        fint = f0*0;
        for i = (1:nelm)
            tmpedof = edof(i, 2:end);
            ed = a(tmpedof);
            ec = [ex(i,:); ey(i,:)];
            [ee, eff] = plan3gs(ec, ed);
            F = F_vect2tens(eff);
            D = mstiff(F,mp);
            es = stresscal(F,mp);
            Ke = plan3ge(ec, t, D, ed, es);
            K(tmpedof, tmpedof) = K(tmpedof, tmpedof) + Ke;
            
            fe = plan3gf(ec, t, ed, es);
            fint(tmpedof) = fint(tmpedof) + fe;
            %disp(max(abs(fe)))
            
        end
        
        G = -fint; % fint-0, fext = 0.
        G(bcn(:,1)) = 0; %vi rï¿½knr med att den blir bra
        %sï¿½tt bc till att vara noll varje varv.
        da = solveq(K, G, bcn);
        a = a + da;
        %fint
        
        G(bcn(:,1)) = 0;
        gplot = [gplot; (norm(G))];
        bcn(right_side, 2) = dbc*0; %för att da inte ska bli för stor
    end
    fplot(n,:) = fint(right_side_nodes);
    uplot(n,:) = a(right_side_nodes);
    if j == maxItr
        disp('reached maxitr')
    end

    an = a;
    %Sn = S;
end
semilogy(gplot)
pause;
plot(uplot, fplot)
pause;
clf
plotpar= [1 4 2];
Ed = extract(edof, a);
eldisp2(ex, ey, Ed, plotpar, 1);
hold on;
plotpar= [1 3 2];
eldraw2(ex, ey, plotpar)