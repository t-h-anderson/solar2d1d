function [ v ] = coeffs2sbp(vin, sim)
% Change from vector of basis coefficients, to piecewise polynomial

v = sim.setup.poly * reshape(vin, sim.setup.pdeg1, sim.setup.nx);

end