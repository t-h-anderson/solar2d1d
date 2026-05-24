function [ out ] = sbptimes(p, q, sim)
%SBPTIMES Takes the product of two sbp functions, returning a new sbp
% function that is the product ignoring the higher order terms.


out = zeros(sim.setup.pdeg1, sim.setup.nx);

for i = 1:sim.setup.nx
    temp = polytimes(p(:,i),q(:,i));
    out(:,i) = polymatch(temp, sim.setup.pdeg1);
end



end

