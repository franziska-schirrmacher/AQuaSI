function [ result ] = shrink( x, thres )
%Computation of the shrink operator for soft thresholding

tmp = abs(x) - thres;
result = (sign(x).* tmp);
setZero = tmp <= 0;
result(setZero) = 0;

end

