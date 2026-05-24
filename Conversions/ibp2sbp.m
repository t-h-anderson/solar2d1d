function [ v ] = sbp2ibp(vin, sim)
% Change from vector of basis coefficients, to piecewise polynomial

v = sim.setup.poly * vin;

end