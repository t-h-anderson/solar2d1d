function [ vout ] = coeffs2space(vin, sim, varargin)
% Change from vector of basis coefficients, to spatial values

v = coeffs2sbp(vin, sim);
if isempty(varargin)
    vout = sbp2space(v, sim);
else
    vout = sbp2space(v, sim, varargin{1});
end

end