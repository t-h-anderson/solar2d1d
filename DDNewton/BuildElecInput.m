
function sim = BuildElecInput(sim)
% Copy input to electrical model
sim.electrical.input = sim.input;
sim.electrical.phys = sim.phys;

% Remove all material data that is only relevant for the optical simulations
f = fields(sim.material);
sim.electrical.material.nsec = sum(sim.material.elecsec);
sim.electrical.material.elecsec = sim.material.elecsec;
for i = 1:length(f)
    if ~strcmp(f{i},'nsec') && ~strcmp(f{i},'elecsec')
        variable = sim.material.(f{i});
        
        if ischar(variable)
            variable  = strtrim(strsplit(variable, '&'));
        end
        
        % Take sections relating to electrical sim
        if length(variable) == sim.material.nsec
            variable = variable(sim.material.elecsec == 1);
        end
        
        if iscell(variable(1))
            variable = strjoin(variable,'&');
        end
        
        sim.electrical.material.(f{i}) = variable;
    end
    
end

% Remove all setup data that is only relevant for the optical simulations
f = fields(sim.electrical.setup);
for i = 1:length(f)
    if ~strcmp(f{i},'nsec') && ~strcmp(f{i},'elecsec')
        variable = sim.electrical.setup.(f{i});
        
        if ischar(variable)
            variable  = strtrim(strsplit(variable, '&'));
        end
        
        % Take sections relating to electrical sim
        if length(variable) == sim.material.nsec
            variable = variable(sim.material.elecsec == 1);
        end
        
        if iscell(variable(1))
            variable = strjoin(variable,'&');
        end
        
        sim.electrical.setup.(f{i}) = variable;
    end
    
end


end


