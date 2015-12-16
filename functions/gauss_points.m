% Lagrar Gausspunkterna i 
% xi_G = [xi_1; ... ; xi_4]
% eta_G = [eta_1; ... ; eta_4]
% samt vikterna, weight = [w1 .. w4]

xi_n = [-1,1,1,-1]';
eta_n = [-1,-1,1,1]';

xi_G = xi_n/sqrt(3);
eta_G = eta_n/sqrt(3);

weight = ones(1,4);
