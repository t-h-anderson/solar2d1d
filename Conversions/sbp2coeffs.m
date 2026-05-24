function [ coeffs ] = sbp2coeffs(sbp, sim)
% Converts a standard basis polynomial to a list of internal basis
% coefficients

coeffs = reshape(sim.setup.ipoly * sbp, sim.setup.np, 1);