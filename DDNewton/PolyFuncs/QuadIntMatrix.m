
function [ integral ] = QuadIntMatrix(func, coeffs, sim)
% Integrate a function * basis poly using Quadrature
      
sim.setup.nx_plot = sim.setup.intdeg;

data = zeros(sim.setup.intdeg * sim.setup.nx, size(coeffs,2));
for i = 1:size(coeffs,2)
    data(:,i) = sim.global.polyvals*coeffs(:,i);
end
f = func(data);

f(isnan(f))=0;
f(isinf(f))=0;

wi = sim.global.wi_quad;

integral = zeros(sim.setup.pdeg, sim.setup.pdeg, sim.setup.nx);
for i = 1:sim.setup.pdeg1
    for j = 1:sim.setup.pdeg1
        integral(i, j,:) = (wi.* sim.global.polyvals2(:,:,i).* sim.global.polyvals2(:,:,j)).'*f;
    end
end

integral = sparse(sim.int.rows, sim.int.cols, reshape(integral,1,[]));

end



