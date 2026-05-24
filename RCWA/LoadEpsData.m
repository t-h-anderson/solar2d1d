



function sim = LoadEpsData(sim)
% Load the permittivity data from data files

% Temporary variables for inputted primary material data
material = strsplit(sim.material.material,'&');
eps1 = zeros(length(sim.setup.nmz), sim.setup.nlambda);

% Temporary variables for inputted secondary material data
eps2 = eps1;
material2 = strsplit(sim.periodic.material2,'&');
periodicregion = 1;

% Loop through each section
for i = 1:sim.material.nsec
    
    % Get the material tag
    secmat = strtrim(material{i});
    
    % If flag is zero, use manually inputted data
    if strcmpi(secmat,'User')
        eps1(sim.setup.matcat == i, :) = sim.material.eps(i);
    % Else load named material
    else
        % Check named material exists
        if exist(secmat, 'file') == 2
            loadedMat = load(secmat);
        elseif exist(strcat(secmat,'.mat'), 'file') == 2
            loadedMat = load(strcat(secmat,'.mat'));
        else
            error('Material Does Not Exist');
        end
        
        % Check permitivity field exists
        if isfield(loadedMat,'eps')
            eps = loadedMat.eps;
            
            % Check type of input
            if isa(eps, 'double') && length(eps) == 1 
                % If a single permittivity then use this for all wavelenghts
                eps1(sim.setup.matcat == i, :) = eps;
            elseif isa(eps, 'double')
                % If an array of permittivities at a range of wavelengths,
                % interpolate for required wavelengths
                eps = interp1(eps(:,1), eps(:,2), sim.setup.nmlambda);
                eps = ones(sum(sim.setup.matcat == i), 1) * eps;
                eps1(sim.setup.matcat == i, :) = eps;
            elseif isa(eps, 'char')
                % If function, assume of form f(nmlambda, VEg)
                eps = str2func(eps);
                try
                    eps = eps(sim.setup.nmlambda, sim.z.Eg(sim.setup.matcat==i));
                    if isnan(eps)
                        error();
                    end
                catch
                    error('Permittivity function not of form eps(nmlambda, VEg)');
                end
                eps1(sim.setup.matcat == i, :) = eps;
            else
                error('Permittivity data not in correct format');
            end
        else
            error('No permittivity data exists for material %s', secmat);
        end
    end
    
    % If section is periodic repeat above procedure but instead look at
    % periodic (secondary) material
    if sim.material.periodicsec(i) == 1
        secmat2 = strtrim(material2{periodicregion});
        
        if strcmp(secmat2,'user') || strcmp(secmat2,'User')
            eps2(sim.setup.matcat == i, :) = sim.periodic.eps2(periodicregion);
        else
            if exist(secmat2, 'file') == 2
                loadedMat = load(secmat2);
            elseif exist(strcat(secmat2,'.mat'), 'file') == 2
                loadedMat = load(strcat(secmat2,'.mat'));
            else
                error('Grating Material Does Not Exist');
            end
            
            if isfield(loadedMat,'eps')
                eps = loadedMat.eps;
                
                if isa(eps, 'double') && length(eps) == 1
                    eps2(sim.setup.matcat == i, :) = eps;
                elseif isa(eps, 'double')
                    eps = interp1(eps(:,1), eps(:,2), sim.setup.nmlambda);
                    eps = ones(sum(sim.setup.matcat == i), 1) * eps;
                    eps2(sim.setup.matcat == i, :) = eps;
                elseif isa(eps, 'char')
                    eps = str2func(eps);
                    
                    eps = eps(sim.phys.nmlambda, sim.z.Eg(sim.setup.matcat==i));
                    eps2(sim.setup.matcat == i, :) = eps;
                else
                    error('Grating permittivity data not in correct format');
                end
            else
                error('No permittivity data exists for grating material');
            end
        end
        
        % Increase counter so next time we look at second grating region
        periodicregion = periodicregion + 1;
        
    end
    
end

% Copy to simulation
sim.zl.eps = eps1;
sim.zl.eps2 =eps2;

end

