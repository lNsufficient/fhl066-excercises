function  d = D_elem(Cinv,a,i,j,k,l)
%Beräknar ett element i materialstyvheten D
% enligt nyhooksk modell
%D_elem(Cinv,a,i,j,k,l)
%Cinv inversen av Cauchytensorn C
%a = [a1 .. a3] a-koefficienter (se manual)
%i,j,k,l index

d=a(1)*Cinv(i,j)*Cinv(k,l)-a(2)*((i==j)*Cinv(k,l)+Cinv(i,j)*(k==l)) ...
    +a(3)*(Cinv(i,k)*Cinv(j,l)+Cinv(i,l)*Cinv(j,k));
