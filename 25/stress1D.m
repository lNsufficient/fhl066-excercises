function S = stress1D(E, ee)
%STRESS1D Summary of this function goes here
%   Detailed explanation goes here

S = E*log(2*ee+1)/(4*sqrt(2*ee+1));

end

