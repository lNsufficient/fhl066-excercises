function coords = coordNbr(ndof)
%COORDNBR Summary of this function goes here
%   Detailed explanation goes here

cols = 3;
coord = [floor(ndof./cols),mod(ndof, cols)]

end

