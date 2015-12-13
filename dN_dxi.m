function val = dN_dxi(i,xi,eta)
%val = dN_dxi(i,xi,eta)

switch i
    case 1
        val = eta-1;
    case 2
        val = 1-eta;
    case 3
        val = 1 + eta;
    case 4
        val = -1-eta;
    otherwise
        disp('Index required: 1 .. 4')
end

val=val/4;
