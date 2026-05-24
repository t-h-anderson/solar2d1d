% Solve the DD equations for a PN junction
% Types of variable
%   - Dimensioned (prepend with units)
% Sorted by type in the sim object

clear();

addpath(genpath('./'));

% Creat the default simulations
sim = DesignJunctionOptDetail();
sim = DesignSimOptDetail(sim);

sim = BuildOptInput(sim);
sim = RCWAinitilisation(sim.optical);
sim = RunOptical(sim);

 if (sim.setup.pol == 0 || sim.setup.pol == 2)
        figure();
        colormap(jet);
        %     % Different graphing options
        contourf(sim.setup.nmx, -sim.setup.nmz, abs(sim.zxl.Ey(:,:,end)), 200,'linestyle','none');
        title('Field: Ey');
        caxis([0,inf])
        axis equal;
        xlim([sim.setup.nmx(1), sim.setup.nmx(end)]);
    end
    
    if (sim.setup.pol == 0 || sim.setup.pol == 1)
        figure();
        colormap(jet);
        %     % Different graphing options
        contourf(sim.setup.nmx, -sim.setup.nmz, abs(sim.zxl.Ep(:,:,end)), 200,'linestyle','none');
        title('Field: Hy');
        caxis([0,inf])
        axis equal;
        xlim([sim.setup.nmx(1), sim.setup.nmx(end)]);
    end




