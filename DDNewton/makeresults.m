function data = makeresults(a, xscale, yscale, sim, nx_plot)

if nx_plot > 0
    x_ref = linspace(-1, 1, nx_plot);
else
    x_ref = -1;
end

% Iterate over each polynomial
data = [];
for i = 1:size(a,2)
    
    % Create grid points dynamically

    x0 = sim.setup.x0(i);
    x1 = sim.setup.x1(i);
       
    if nx_plot > 0
        x_t = linspace(x0, x1, nx_plot);
        sigvals = polyval(a(:,i), x_ref);
    else
        x_t = x1;
        sigvals = a(i);
    end
            

   data = [data,[xscale .* x_t;(yscale * sigvals)]];
    
    x0 = x1;
end

data = data.';

end

