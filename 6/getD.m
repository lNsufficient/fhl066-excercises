function D = getD(E, v)
%GETD Summary of this function goes here
%   Detailed explanation goes here
D = E/(1-v*v)*[1, v, 0; v, 1, 0; 0, 0, (1-v)/2];

end

