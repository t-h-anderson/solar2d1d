function [ fp ] = func2sbp(f, sim, varargin)
% func2poly(func, pdeg1, h)
% Turns a function 'func' evaluated over (0,sum(h))
% into a piecewise poly in sbp, i.e. each poly defined
% on the standerd basis element

if ~isempty(varargin)
    mesh = sim.setup.mesh(cell2mat(varargin)==1);
    nx = length(mesh);
else
    mesh = sim.setup.mesh;
    nx = sim.setup.nx;
end

pdeg1 = sim.setup.pdeg1;
h = mesh/sum(mesh); % Scale between 0 and 1

% Get locations for polynomial eval on element
xbar = fliplr(lgnodes(pdeg1-1));

fp = zeros(pdeg1, nx);
x0 = 0;
for i = 1:nx
    x1 = x0 + h(i);
    
    % Find x locations between x0 and x1
    x = x0+(xbar+1)*(x1-x0)/2;
    fvals = f(x);
    if length(fvals) == 1
        fvals = fvals * ones(1,pdeg1);
    end
    
    % Fit the polynomial
    fp(:,i) = polyfit_nowarn(xbar, fvals, pdeg1-1);
    
    %Move to next element
    x0=x1;
end

end

