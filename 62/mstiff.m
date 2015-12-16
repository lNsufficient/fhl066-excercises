function D = mstiff(F,mp)
%beräknar materialstyvheten D
%D = mstiff(F,mp)
% mp = [E,nu]
E=mp(1);    nu=mp(2);

J = det(F);
C=F'*F;
Cinv=inv(C);
K=E/(3*(1-2*nu));
G=E/(2*(1+nu));

a=zeros(3,1);
a(1)=K*J^2+(2*G/9)*(J^(-2/3))*trace(C);
a(2)=(2*G/3)*J^(-2/3);
a(3)=(G/3)*trace(C)*J^(-2/3)-(K/2)*(J^2-1);

D1111=D_elem(Cinv,a,1,1,1,1);
D2222=D_elem(Cinv,a,2,2,2,2);
D1212=D_elem(Cinv,a,1,2,1,2);
D1122=D_elem(Cinv,a,1,1,2,2);
D1112=D_elem(Cinv,a,1,1,1,2);
D2212=D_elem(Cinv,a,2,2,1,2);

D = [D1111/2, D1122,    D1112;
    0       , D2222/2,  D2212;
    0   ,       0   ,   D1212/2];

D=D+D';
