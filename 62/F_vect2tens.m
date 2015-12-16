function F_tens = F_vect2tens(F_vect)

% F_vect deformationsgradient - vektor

F_tens=zeros(3);

F_tens(1,1) = F_vect(1); %dx/dx0
F_tens(2,2) = F_vect(4); %dy/dy0

F_tens(2,1) = F_vect(3); %dy/dx0
F_tens(1,2) = F_vect(2); %dx/dy0

% Plan töjning
F_tens(3,3) = 1; %dz/dz0
% Resten 0
