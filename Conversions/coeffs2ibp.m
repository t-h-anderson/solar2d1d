function [ v ] = coeffs2ibp(vin, sim)
% Change from vector of basis coefficients, to piecewise polynomial

v = reshape(vin, sim.setup.pdeg1, sim.setup.nx);

end