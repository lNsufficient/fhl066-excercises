function val = dN_deta(i,xi,eta)
%val = dN_deta(i,xi,eta)

switch i
    case 1
        val = xi-1;
    case 2
        val = -1-xi;
    case 3
        val = 1 + xi;
    case 4
        val = 1-xi;
    otherwise
        disp('Index required: 1 .. 4')
end

val=val/4;
