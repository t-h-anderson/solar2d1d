function [ integral ] = QuadIntVect(func, coeffs, sim)
% Integrate a function * basis poly using Quadrature
      
%sim.setup.nx_plot = sim.setup.intdeg;

data = zeros(sim.setup.intdeg * sim.setup.nx, size(coeffs,2));
for i = 1:size(coeffs,2)
    data(:,i) = sim.global.polyvals*coeffs(:,i);
end
f = func(data);

f(isnan(f))=0;
f(isinf(f))=0;

wi = sim.global.wi_quad;

integral = zeros(sim.setup.nx, sim.setup.pdeg1);
for i = 1:sim.setup.pdeg1
    integral(:,i) = (wi.* sim.global.polyvals2(:,:,i)).'*f;
end

integral = reshape(integral.',[],1);

end


