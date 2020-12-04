function [ matrix ] = vectorToMatrix( vector,R,C,D)
%Function converts a vector into a matrix. The storage is made in a
%column-wise manner

if nargin < 4
    matrix = reshape(vector,R,C);
else
    matrix = reshape(vector,R,C,D);
end

end

