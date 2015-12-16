function S = stresscal(F,mp)
%ber√§knar 2 PK-sp√§nningen S enligt nyhooksk modell
% F ‰r def.gradient pÂ matrisform
% S = [S_xx, S_yy, S_xy]' = stresscal(F,mp)
% mp = [E,nu]
E=mp(1);    nu=mp(2);

J = det(F);
C=F'*F;
Cinv=inv(C);

%G = 3*K*(1-2*nu)/(2*(1+nu));
K=E/(3*(1-2*nu));
G=E/(2*(1+nu));

S_tens = (K/2)*(J^2-1)*Cinv+G*(J^(-2/3))*(eye(3)-(trace(C)*Cinv/3));
S = [S_tens(1,1);S_tens(2,2);S_tens(1,2)];
