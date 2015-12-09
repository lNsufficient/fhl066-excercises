function Ke=bar3ge_log(ec,ep,ed,ee)

%Ke=bar3ge(ec,ep,ed,es,eps)
% Ke elementstyvhet
% ---
% ec odeformerade koordinater [x1 x2; y1 y2 ; z1 z2]
% ep materialparametrar [E A0]
% ed f�rskjutningar [a1 a2 ... a6]
% es normalkraft
% ee greent�jning

x0 = bar(ec);
%l0=sqrt(x0'*x0);
l_0=norm(x0);
E = ep(1); A_0=ep(2);
x = bar_def(ec,ed);
dedu = [x',-x']/(l_0^2);
dxdu = [eye(3), -eye(3); -eye(3), eye(3)];

dsde = dmat1D(E, ee);
S = stress1D(E, ee);

Ke = (A_0/l_0)*(dsde*[x;-x]*dedu+S*dxdu);

% 
% 
% 
% 
% u=(ed(4:6)-ed(1:3)); %ber�knar f�rskjutning
% 
% E=ep(1); A0 = ep(2);
% 
% %Enligt 2.35-2.37:
% M=x0*x0';
% K0 = E*A0/l0^3*[M -M; -M M];  
% 
% U=x0*u'+u*x0'+u*u';
% 
% Ku = E*A0/l0^3*[U -U; -U U];
% 
% Ksigma = es/l0*[eye(3) -eye(3); -eye(3) eye(3)];
% 
% Ke = K0 + Ku + Ksigma; %Enligt 2.34