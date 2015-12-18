clear;
perturb_switch = 0;
%0. ingen st?rn
%2. st?rn i lasten

data_e8 %denna kör automatiskt data.m

plot_a = [];
plot_f = [];

lambda=0;
maxItr = 100;

if perturb_switch == 2
    P(top_dof) = P_end*cosd(12);
    P(top_dof-1) = P_end*sind(12);
end
G = P*0;
f_int=P*0;
a = P*0;
K = zeros(ndof, ndof);

%TOL=norm(P_end/nbr_steps)*1e-2;
SCALE=1;
l=24e0;
l_0=l;
TOL = 1e-2;

n_end=400;

psy=0; %given in task

n=0;

top_dof=find(P);

old_a = a; 

USE_HOOKE = 0;
continuous_plot = 1;

plot_res=[];

while n < n_end
    n=n+1;
    a_i = a;
    lambda_i = lambda;
 
    res = TOL+1;
    i = 1;
    G=G*0;
    while (res > TOL) %This could cause a problem, res = 0 first itr?
        K=K*0;
        f_int=f_int*0;

        for j = 1:nelm
            index_dof=edof(j,2:end); %de frihg stången gränsar till
            ec=[ex(j,:); ey(j,:)];
            ed=a_i(index_dof);
            [ee,eff]= plan3gs(ec,ed);
            
            F = F_vect2tens(eff);
            D = mstiff(F,mp);
            es = stresscal(F,mp);
            
            Ke = plan3ge(ec,t,D,ed,es);
            K(index_dof,index_dof)=K(index_dof,index_dof)+Ke;
            
            %ef= plan3gf(ec, t, ed, es);
            %f_int(index_dof) = f_int(index_dof) + ef;
        end
        
        da_G = solveq(K,-G,bc);
        da_P = solveq(K,P,bc);
        
        delta_a=a_i-a;
        delta_lambda=lambda_i-lambda; 
        [a1, a2, a3, a4, a5] = calc_a(P, delta_a, delta_lambda, da_P, da_G, l, psy);
        
        if i > 1
            p=a2/a1;        q=a3/a1;
            d_lambda_test=-p/2+sqrt((p/2)^2-q)*[1; -1];
            if (a4+a5*d_lambda_test(1) > a4+a5*d_lambda_test(2))
                d_lambda_test=d_lambda_test(1);
            else
                d_lambda_test=d_lambda_test(2);
            end
        else
            delta_a_n = a - old_a;
            s=sign(delta_a_n'*da_P);
            if (n == 1)
                s = 1;
            end
            d_lambda_test=s*sqrt(-a3/a1);
            s=1;
        end   
        if ~isreal(d_lambda_test)
            keyboard
            a_i = old_a;
            
            l = l/2;
            continue;
        else
            d_lambda = d_lambda_test;
            l = min(l*2,l_0);
        end
            
        lambda_i = lambda_i + d_lambda; 
        
        a_i = a_i + da_G + da_P*d_lambda;
 
        for j = 1:nelm
            index_dof=edof(j,2:end); %de frihg stången gränsar till
            ec=[ex(j,:); ey(j,:)];
            ed=a_i(index_dof);
            [ee,eff]= plan3gs(ec,ed);
            
            F = F_vect2tens(eff);
            D = mstiff(F,mp);
            es = stresscal(F,mp);
            
            %Ke = plan3ge(ec,t,D,ed,es);
            %K(index_dof,index_dof)=K(index_dof,index_dof)+Ke;
            
            ef= plan3gf(ec, t, ed, es);
            f_int(index_dof) = f_int(index_dof) + ef;
        end
        
        G=f_int-lambda_i*P;
        G(bc(:,1)) = 0;
        res=norm(G);
        plot_res = [plot_res, res];
        i=i+1;
        
        
        if i == maxItr
        %    keyboard;
            break;
        end
        
        
    end
    
    
    old_a = a;
    a=a_i;

    if n==30
    oldLambda = lambda-lambda_i;
    end
    
    
    if (n > 30 &&((lambda_i-lambda)/oldLambda > 2))
        %keyboard
    end

    lambda=lambda_i;
   
    plot_f(n) = lambda*P_end;
    plot_a(n) = a(top_dof);
    
    if mod(n,500)==0
        plot(plot_a, plot_f)
        disp('========')
    end
    
    if (continuous_plot && ~mod(n, 100))
        clf;
        plotpar= [1 4 3]; 
        Ed = extract(edof, a);
        subplot(2,1,1);
        eldisp2(ex, ey, Ed, plotpar, 1);
        subplot(2,1,2);
        plot(plot_a, plot_f)
        pause;
    end
    
    [n,i-1]
    
end


%%
%plot(abs(plot_a), abs(plot_f))
d_lambda = (plot_f(2:end) - plot_f(1:end-1))/P_end;
semilogy(plot_res)
plot(real(d_lambda))
plot(plot_a, plot_f)
xlabel('f�rskjutning / meter')
ylabel('kraft / Newton')
A = zeros(length(plot_a), 3);
F = zeros(length(plot_f), 3);
A(:,perturb_switch+1) = plot_a;
F(:,perturb_switch+1) = plot_f;
hold off;
clf;
%%
clf;
plotpar= [1 4 3]; 
Ed = extract(edof, a);
eldisp2(ex, ey, Ed, plotpar, 1);
hold on;
plotpar= [1 2 2];
eldraw2(ex, ey, plotpar)