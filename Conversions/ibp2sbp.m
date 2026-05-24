function [ v ] = ibp2sbp(vin, sim)
% Convert an individual-basis-polynomial representation to the standard
% basis polynomial form by left-multiplying by sim.setup.poly.

v = sim.setup.poly * vin;

end