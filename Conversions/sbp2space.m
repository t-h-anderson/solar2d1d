function [data] = sbp2space(p, sim, varargin)
% Convert sbp to spacial values



if isempty(varargin)
    nx_plot = sim.setup.nx_plot;
    x_ref = linspace(-1, 1, nx_plot);
else
    x_ref = sim.local.xloc_qdatauad';
    nx_plot = length(x_ref);
end

data = zeros(size(p,2) * nx_plot, 2);

for i=1:size(p,2)
    if isfield(sim.setup,'x_l')
        x = linspace(sim.setup.x_l(i), sim.setup.x_r(i), nx_plot);
    else
        x = linspace(sim.setup.x0(i), sim.setup.x1(i), nx_plot); 
    end
    y = polyval(p(:, i), x_ref);
    
    data((i - 1)*nx_plot+1:i*nx_plot,:) = [x; y]';
end

end