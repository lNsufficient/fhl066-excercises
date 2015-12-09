function D = dmat1D(E, ee)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

D = E/4*(2*1/(2*ee+1)*1/(sqrt(2*ee+1))+2*log(2*ee+1)*(-1/2)*(2*ee+1)^(-3/2));
end

