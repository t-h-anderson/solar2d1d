function makeplot(a, plottitle, xscale, yscale, sim, plotno, clearplot, colour, nx_plot)

if nx_plot > 0
    x_ref = linspace(-1, 1, nx_plot);
else
    x_ref = -1;
end

% Set figure
subplot(sim.setup.plotlayout(1), sim.setup.plotlayout(2), plotno);

% Clear plot if requested
hold on
if clearplot
    hold off
    plot(0,0)
end


% Iterate over each polynomial

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
            
    if sim.setup.lp       
        plot(xscale .* x_t(sigvals > 0), log10(yscale*sigvals(sigvals > 0)), colour);
        if(colour(end)~='-')
            plot(xscale .* x_t(sigvals < 0),log10(-yscale*sigvals(sigvals < 0)), strcat(colour,':'));
        end
        hold on
    else
        plot(xscale .* x_t,(yscale * sigvals),colour);
        hold on
    end
    
    x0 = x1;
end

if ~strcmp(plottitle, '')
    title(plottitle);
end

xlim(xscale*[sim.setup.x0(1),sim.setup.x1(end)])

hold off



end

