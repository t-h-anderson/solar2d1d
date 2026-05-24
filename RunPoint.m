% Solve the DD equations for a PN junction
% Types of variable
%   - Dimensioned (prepend with units)
% Sorted by type in the sim object

pointno = 1744;

addpath(genpath('./'));

% Creat the default simulations
%sim = DesignJunctionOpticTest();
%sim = DesignSimOpticTest(sim);

for i = 1:size(optimResult.paramDefCell,1)
    fieldname = optimResult.paramDefCell{i,1};
    params.(genvarname(fieldname)) = optimResult.allTestedMembers(i,pointno);
end

sim = DesignJunction();
sim = DesignSim(sim);

% Put in mapping here
sim.material.nmSec(4) = params.nmLi;
sim.material.nmSec(6) = params.nmLg;
sim.optical.periodic.zeta = params.zeta;
sim.input.nmLx = params.nmLx;

switch 2
    case 1
        Eg0 = params.Eg0;
        Eg1 = params.Eg1;
        
        m = num2str((Eg1-Eg0)/(params.nmLi+40));
        f1 = strcat('@(x)', m, '*x*20');
        f2 = strcat('@(x) ',num2str(m),'*(',num2str(params.nmLi),'*x + 20)');
        f3 = strcat('@(x) ',m,'*(20*x +', num2str(params.nmLi), '+ 20)');
        
        sim.material.VEg1 = strcat('0&0&',f1,'&',f2,'&',f3,'&0&0');
        sim.material.VEg0 = Eg0;
    case 2
        Eg_p = params.Eg_p;
        Eg_i = params.Eg_i;
        Eg_n = params.Eg_n;
        
        sim.material.VEg0 = Eg_i;
        sim.material.VEg1 = strcat('0&0&',num2str(Eg_p-Eg_i),'&0&',num2str(Eg_n-Eg_i),'&0&0');
        
    case 3
        %sim.material.VEg1 = strcat('0&0&0&@(x)',num2str(A),'*(0.5*(1+sin(',num2str(kappa),'*2*pi*x - ',num2str(phi),'*2*pi))).^',num2str(alpha),'&0&0&0');
end
%sim.material.icm3ND = [0,0,10^params.Nd,0,-10^params.Na,0,0];

sim = BuildOptInput(sim);

% Scale mesh to keep same resolution
%sim.optical.setup.nslices(3) = ceil(params.nmLp/2);
sim.optical.setup.nslices(4) = ceil(params.nmLi/2);
%sim.optical.setup.nslices(5) = ceil(params.nmLn/2);

nmdz = sim.electrical.setup.nmdz; %
if nmdz ~= 0
    sim.electrical.setup.nx_sec = ceil(sim.material.nmSec/nmdz);
end


if sim.optical.setup.optical_toggle
    % Copy input to optical
    sim.optical = RCWAinitilisation(sim.optical);
    sim.optical = RunOptical(sim.optical);
    sim.input.icm3isG = sim.optical.results.icm3isG;
end

if sim.electrical.setup.electrical_toggle
    % Run electrical model
    sim = BuildElecInput(sim);
    sim.electrical = RunElectrical(sim.electrical);
    
    if sim.electrical.setup.plotsolve > 0
        figure(2)
        hold off
        plot(abs(sim.electrical.results.V), sim.electrical.results.J, '.')
        hold on
        plot(abs(sim.electrical.results.V), sim.electrical.results.P, '.')
        ylim([0, 1.1* max(max([sim.electrical.results.P, sim.electrical.results.J]))])
    end
    
    % disp(sim.electrical.results)
    disp(sim.electrical.results.Wim2Pmax/10);
    
end



