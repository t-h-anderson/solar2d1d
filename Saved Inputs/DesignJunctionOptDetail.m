function [sim] = DesignJunctionOptDetail()
%JUNCTION Create a junction object. If specified property is length 1,
%then it will be made the correct length automatically.

sim.input = struct(...
    'nmLx', 500, ...       % Grating period
    'VVext', 0.0, ...       % External voltage
    'icm3isG', 1e20, ...   % Generation Rate (value, function, file, or -inf for results from RCWA)
    ...% Junction Material Specification (used if materlial = 1)
    'Eg0', 1.6, ...         % Baseline Bandgap
    'A', 0., ...            % Perturbation Amplitude
    'kappa', 3702, ...         % Periods
    'phi', 0.75, ...        % Phase
    'alpha', 5 ... % Shaping Parameter
    );

sim.material = struct(...
    'nsec', 7, ...         % Number of independent section in study
    'elecsec', [0, 0, 1, 1, 1, 0, 0], ... % Flags for which sections to include in the electrical simulation
    'periodicsec', [0, 0, 0, 0, 0, 1, 0], ... % Flags for which sections to include in the electrical simulation
    'nmSec', [200, 75, 20, 200, 20, 50, 100], ... % Lengths of the sections
    'icm3ND', [0, 0, 1e16, 0, -1e16, 0, 0], ... % Doping of the sections
    ... % Material properties will be overwritten if 'material' not set to 0. If material has a range of possible bandgaps, will tailor it to specified Eg.
    'material', 'User & User & User & User & User & User & User' , ... % Loads a material from the database, 0 for user input %    'material', 'Air & AZO & CIGS & CIGS & CIGS & AZO & Silver' , ... % Loads a material from the database, 0 for user input
    'VEg0', 1.3, ... % Baseline bandgap
    'VEg1', 0,...   % Bandgap perturbation
    ... % Electrical Properties
    'icm3Nc', 2.5e20, ... % Conduction band DOS
    'icm3Nv', 2.5e20, ... % Valence band DOS
    'cm2iVismun', 15, ...  % Electron mobility
    'cm2iVismup',  15, ...  % Hole mobility
    'VChi',  4, ...% Electron Affinity
    'epsdc', 11.68, ... % DC relative permittivity of bulk silicon
    ... % Recombination Parameters
    'icm3isalpha', 1e-12, ... % Radiative recombination constant 1/(cm^3 s), related to beta via beta = alpha/ni^2
    'istausrhn', 1e-9, ... % SRH recombination constant 1/s
    'istausrhp', 1e-9, ... % SRH recombination constant 1/s
    'Cn', 1e-28, ... % Auger recombination constant
    'Cp', 1e-28, ... % Auger recombination constant
    'ET', 0.5, ... % Relative position of trap level in bandgap
    ... % Optical Properties
    'eps', [1+1e-9i, 3.33+0.016i, 9.5+1.25i, 9.5+1.25i,9.5+1.25i, 3.33+0.016i, -22+0.4i] ... % Optical Permittivity
    );

sim.optical.periodic = struct(...
    ... % Takes a vector of inputs. One entry per periodic region specified above
    'zr_flag', [1], ... Flag to toggle between 1 zeta, and 2 relief
    'zeta', '@(z) 0.5', ... %'@(z) (1/pi) * acos(z/80)', ...         % Graing Duty Cycle
    'relief', '@(zeta) cos(2*pi*zeta)', ... %%& @(zeta) 1-cos(2*pi*zeta)', ...
    ...&'type', 1, ...                                  % Grating Shape (0 for rect., 1 for sine, 2 for triangle/trap) or function (to do)
    'material2', 'User', ...                      % Grating material ( 1 for silver, const otherwise)
    'eps2', [-22+0.4i] ...                                   % Grating metal permittivity
    );

% Constants
sim.phys.Cq = 1.6021766e-19;
sim.phys.Jshbar = 1.0545718e-34;
sim.phys.VVth = 300 * 8.61733034e-5;
sim.phys.CiVimeps0 = 8.85418782e-12; % converted im to icm
sim.phys.CiVicmeps0 = sim.phys.CiVimeps0 / 100; % converted im to icm

% Units needed here
sim.phys.mkgis2iA2mu0=4*pi*10^(-7);           % permeability of free space
sim.phys.Oeta0 = 376.73031346177; % Impedence of free space
sim.phys.misc = 299792458;                 % speed of light
sim.phys.m2kgish = 6.6260696e-34;       % Planck constant
sim.phys.Enormconst = sqrt(1/(sim.phys.CiVimeps0 * sim.phys.misc));% sqrt(2*sim.phys.Oeta0); % Scale E by this and H by 1/this such that incident power density is 1 W/m^2

% Scaling (These can be changed)
sim.phys.nmLs = sim.material.nmSec * sim.material.elecsec.'; % in nm
sim.phys.cmLs = cm_from_nm(sim.phys.nmLs); % in nm
sim.phys.icm3Ns = max([max(abs(sim.material.icm3ND)),1e14]); % in cm^-3
sim.phys.VVs = sim.phys.VVth; % in V
sim.phys.cm2iVismus = 1.5; % in cm^2V^-1s^-1
sim.phys.icm3isGs = sim.phys.cm2iVismus * sim.phys.VVs * sim.phys.icm3Ns/sim.phys.cmLs^2;

sim.phys.cm3isalphas = (sim.phys.cm2iVismus * sim.phys.VVs)/(sim.phys.cmLs^2 * sim.phys.icm3Ns);
sim.phys.isTaus = sim.phys.cmLs^2/(sim.phys.cm2iVismus * sim.phys.VVs);
sim.phys.cm6isCs = (sim.phys.cm2iVismus * sim.phys.VVs)/(sim.phys.cmLs^2 * sim.phys.icm3Ns^2);

%%%%%%%%%%%%%%%%%%%%%%
% Do not edit below this box %
%%%%%%%%%%%%%%%%%%%%%%
nsec = sim.material.nsec;
if length(sim.material.nmSec) ~= nsec
    error('Lengths not defined for each section of junction');
end

% Convert single inputs to a constant input for each section
sim.material = convert12nsec(sim.material, sim.material.nsec);
sim.optical.periodic = convert12nsec(sim.optical.periodic, sum(sim.material.periodicsec));
sim.optical.periodic.nsec = sum(sim.material.periodicsec);

% Other interal scaling factors
sim.phys.cmLs = cm_from_nm(sim.phys.nmLs); % Ls in cm
sim.phys.cm2isDs = sim.phys.cm2iVismus * sim.phys.VVs; %in cm^2s^-1
sim.phys.mAicm2Js = 1000 * sim.phys.Cq * sim.phys.VVs * sim.phys.cm2iVismus * sim.phys.icm3Ns/sim.phys.cmLs; % in mA cm^-2
sim.phys.icm2isGRs = sim.phys.VVs * sim.phys.cm2iVismus * sim.phys.icm3Ns/sim.phys.cmLs^2; % in cm^-3s^-1
sim.phys.lambda2s = sim.phys.CiVicmeps0 * sim.phys.VVth/(sim.phys.cmLs^2 * sim.phys.Cq * sim.phys.icm3Ns); % Dimensionless (needs multiplying by eps)
sim.phys.Qscale = sim.phys.Oeta0 * sim.phys.CiVimeps0/(sim.phys.Jshbar * sim.phys.Enormconst^2);

end

% Convert a length 1 input to an nsec length
function material = convert12nsec(material, nsec)

f = fields(material);
% Loop through the fields
for i = 1:length(f)
    % Get the material property
    prop = material.(f{i});
    if ischar(prop) % if we have a character string split at '&'
        propcells = strsplit(prop,'&');
        if length(propcells) == 1 % If length is 1, then duplicate to nsec
            prop = strtrim(propcells{1});
            p2 = prop;
            for j = 1:nsec-1
                prop = strcat(prop, '& ', p2);
            end
            material.(f{i}) = prop;
        elseif length(propcells) ~= nsec
            error('Input "%s" incorrect length', f{i});
        else
            prop = strtrim(propcells{1});
            if ~strcmp(prop(1),'@') && ~isempty(str2num(prop(1)))
                 prop = strcat('@(x)',prop);
            end
            
            for j = 2:nsec
                p2 = strtrim(propcells{j});
                if ~strcmp(p2(1),'@') && ~isempty(str2num(p2(1)))
                    p2 = strcat('@(x)',p2);
                end
                prop = strcat(prop, '& ', p2);
            end
            material.(f{i}) = prop;
        end
    else 
        if length(prop) == 1 && ~strcmp(f{i}, 'nsec') % If we have one number, duplicate to nsec (ignore nsec)
            prop = prop*ones(1,nsec);
            material.(f{i}) = prop;
        elseif length(prop) ~= nsec && ~strcmp(f{i}, 'nsec')
            error('Input "%s" incorrect length', f{i});
        end
    end
end

end

