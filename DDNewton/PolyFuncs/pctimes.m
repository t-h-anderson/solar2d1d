function [ M ] = pctimes(c, v)
% Function which multiplies the vector of poly coefficients by piecewise constant function.

pdeg = size(v,1)/length(c);
C = diag(repmat(c,1,pdeg));
M = C*v;

end