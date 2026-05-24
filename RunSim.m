clear();
clc;
addpath(genpath('./'));

% Create the simulations
sim = DesignJunction();
sim = DesignSim(sim);


if sim.optical.setup.optical_toggle
    % Copy input to optical model
    sim = BuildOptInput(sim);
    
    % Initialise optical model 
    sim.optical = RCWAinitilisation(sim.optical);
    
    % Run optical model
    sim.optical = RunOptical(sim.optical);
       
    % Transfer results to electrical model
    sim.input.icm3isG = sim.optical.results.icm3isG;
end
disp('Optical done');
if sim.electrical.setup.electrical_toggle
    
    % Build electrical model
    sim = BuildElecInput(sim);
    
    % Run electrical model
    sim.electrical = RunElectrical(sim.electrical);
   
end
disp('Electrical done');


