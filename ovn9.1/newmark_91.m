clear;
data_e9;

f_int = zeros(ndof,1);  

%i det h‰r fallet ‰r g bara den interna kraften, ty C = 0.
for j = 1:nelm
            index_dof=edof(j,2:end); %de frihg st√•ngen gr√§nsar till
            ec=[ex(j,:); ey(j,:)];
            ed=a(index_dof);
            [ee,eff]= plan3gs(ec,ed);
            
            F = F_vect2tens(eff);
            D = mstiff(F,mp);
            es = stresscal(F,mp);
            
            %Ke = plan3ge(ec,t,D,ed,es);
            %K(index_dof,index_dof)=K(index_dof,index_dof)+Ke;
            
            ef= plan3gf(ec, t, ed, es);
            f_int(index_dof) = f_int(index_dof) + ef;
end   

%1 - initial conditions
upp_0 = Minv*(f_0-g(u_0, up_0));

n = 0;
n_end = 40;
while(n < n_end)
    n = n+1;
    %2 - prediction step
    upp_np1 = upp_n;
    up_np1 = up_n + h*upp_n;
    u_np1 = u_n + h*up_n+1/2*h^2*upp_np1;
    
    res = eps_r + 2;
    du_norm = eps_u + 3;
    while (res > eps_r || du_norm > eps_u)
        %3 - residual calculation
        r = f_np1 - M*upp_np1-g(u_np1,up_np1);
        
        %4 - system matrices and increment correction
        K = getK();
        C = K*0;  %vi har antagit att det inte beror pÂ hastigheter
        K_star = K + gamma/(beta*h)*C + 1/(beta*h^2)*M;
        du = solveq(K_star, r, bc);
        u_np1 = u_np1 + du;
        up_np1 = up_np1 + gamma/(beta*h)*du;
        upp_np1 = upp_np1 + 1/(beta*h^2)*du;
        res = norm(r);
        du_norm = norm(du);
    end
end