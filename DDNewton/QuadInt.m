function [ integral ] = QuadInt(func, coeffs, sim)
% Integrate a function * basis poly using Quadrature
      
sim.setup.nx_plot = sim.setup.intdeg;

px  = sim.local.polyvals;

data = zeros(sim.setup.intdeg * sim.setup.nx, size(coeffs,2));
for i = 1:size(coeffs,2)
    data(:,i) = sim.global.polyvals*coeffs(:,i);
end
f = func(data);

wi = sim.global.wi_quad;


%integral = wi*f;

integral = zeros(sim.setup.nx, sim.setup.pdeg);
for i = 1:sim.setup.pdeg1
    
    integral(:,i) = (wi.* sim.global.polyvals2(:,:,i)).'*f;
    
    %integral(index1:index2) = sim.setup.poly*integral(index1:index2);
end

integral = reshape(integral.',[],1);

end