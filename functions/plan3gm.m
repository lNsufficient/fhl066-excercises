function Me = plan3gm(ec,t,rho)
ex=ec(1,:);
ey=ec(2,:);
A = triarea(ex,ey);
er=[0,1,0];
es=[0,0,1];
N1 = @(x,y) (er(2)*es(3)-er(3)*es(2)+(es(2)-es(3))*x+(er(3)-er(2))*y);%/(2*A);
N2 = @(x,y) (er(3)*es(1)-er(1)*es(3)+(es(3)-es(1))*x+(er(1)-er(3))*y);%/(2*A);
N3 = @(x,y) (er(1)*es(2)-er(2)*es(1)+(es(1)-es(2))*x+(er(2)-er(1))*y);%/(2*A);
%kesboard
%for i=1:3
    %[xy,w] = dunavant_rule(2,3);
    [xy,w] = dunavant_rule(2);
%end
Me=zeros(6);
for i=1:3
    
            f1=N1(xy(1,i),xy(2,i));
            f2=N2(xy(1,i),xy(2,i));
            f3=N3(xy(1,i),xy(2,i));
    N=[f1,0,f2,0,f3,0;
        0,f1,0,f2,0,f3];
    Me=Me+N'*N*w(i);
end

Me=Me*A*t*rho;