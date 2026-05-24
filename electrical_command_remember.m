clear all
close all;
clc;


clear all
close all
clc;
addpath(genpath('./'))

% get default DE parameters
DEParams = DesignDE();
Nx=100;
nsec=9;
nslices=[1,1,1,1,25,100,1,50,1];
nmSec=[100,110,100, 70, 50, 200,10, 150, 100];
elecsec= [0,0,0,0,1,1,0,0,0];
periodicsec= [0,0,0,0,0,0,0,1,0];
zr_flag=1;
zeta=0.5;
relief=@(zeta) cos(2*pi*zeta);
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
sim.phys.Enormconst = sqrt(1/(sim.phys.CiVimeps0 * sim.phys.misc)); % sqrt(2*sim.phys.Oeta0); % Scale E by this and H by 1/this such that incident power density is 1 W/m^2

% Scaling (These can be changed)
sim.phys.nmLs = nmSec * elecsec.'; % in nm












nslices=nslices.*abs(sign(nmSec));
Nz=sum(nslices);
nmz=zeros(Nz,1);
nmdz=zeros(Nz,1);
matcat=zeros(Nz,1);
z0=0;
z1=0;
index1=1;

elecmask=zeros(Nz,Nx);
nmdzmask=zeros(Nz,Nx);
for sec=1:length(nslices)
    index2 =index1+nslices(sec)-1;
    nmdz(index1:index2) =nmSec(sec)/nslices(sec);
    z1=z1+nmSec(sec);
   

if z1 ~= z0
        newnmz = linspace(z0 + nmdz(index1)/2, z1 - nmdz(index1)/2, nslices(sec));
    else
        newnmz = [];
    end
    nmz(index1:index2) = newnmz;
    
    matcat(index1:index2) = sec;
    
    elecmask(index1:index2, :) = elecsec(sec);
    nmdzmask(index1:index2, :) = nmSec(sec)/nslices(sec);
   
    z0=z1;
    index1 = index2+1;
end   
    sec = 1;
sim.setup.nmz0 = 0;
while sec <= length(elecsec) &&elecsec(sec)==0
    sim.setup.nmz0 = sim.setup.nmz0 + nmSec(sec);
    sec = sec + 1;
   
end

sim.setup.nmz1 = sim.setup.nmz0;
while sec <= length(elecsec) && elecsec(sec)==1
    sim.setup.nmz1 = sim.setup.nmz1 + nmSec(sec);
    sec = sec + 1;
    
end

sim.setup.nmLz = sim.setup.nmz1 - sim.setup.nmz0;
sim.setup.Lz = sim.setup.nmLz / sim.phys.nmLs;
disp(sim.setup.nmLz);
sim.setup.nmz = nmz;
sim.setup.nmdz = nmdz;
sim.setup.Nz = Nz;
sim.setup.matcat = matcat;

sim.zx.elecmask = elecmask;
sim.zx.nmdzmask = nmdzmask;

sim.z.elecmask = elecmask(:,1);
sim.z.nmdzmask = nmdzmask(:,1);
 
zetas = linspace(-0.5,0.5,Nx);
zetalist = [];
periodiccount =1;
loc = 0;
maxlength = 2;

for sec = 1:nsec
    
    % If periodic section
    if periodicsec(sec) == 1
        
        % If def by zeta
        if zr_flag(periodiccount) == 1
            
            % If number
            if isa(zeta(periodiccount), 'double')
                % Set all layers to number
                for i = 1:nslices(sec)
                    zetas = [-0.5*zeta(periodiccount), 0.5*zeta(periodiccount)];
                    if(zetas(1) ~= -0.5)
                        zetas = [-0.5, zetas];
                    end
                    if(zetas(end) ~= 0.5)
                        zetas = [zetas, 0.5];
                    end
                    zetalist{i+loc} = zetas';
                    maxlength = max(maxlength, length(zetas));
                end
            else
                zetacells = strsplit(zeta,'&');
                zetaFn = zetacells{periodiccount};
                
                
                
                zetaFn = str2func(zetaFn);
                
                dz = 1/ nslices(sec);
                z=1-linspace(0+0.5*dz, 1-0.5*dz,nslices(sec));
                for i = 1:nslices(sec)
                    zetas = [-0.5*zetaFn(z(i)), 0.5*zetaFn(z(i))];
                    if(zetas(1) ~= -0.5)
                        zetas = [-0.5, zetas];
                    end
                    if(zetas(end) ~= 0.5)
                        zetas = [zetas, 0.5];
                    end
                    zetalist{i+loc} = zetas';
                    maxlength = max(maxlength, length(zetas));
                end
            end
        elseif zr_flag(periodiccount) == 2 % Specify sim.z.zeta via input surface
            
            % Load function for section
            reliefcells = strsplit(relief,'&');
            reliefFn = reliefcells{periodiccount};
            reliefFn = str2func(reliefFn);
            zeta = linspace(-0.5, 0.5, Nx);
            
            % Surface g(x)
            relief = reliefFn(zeta);
            
            % Normalized z positions in periodic region
            dz = 1/ nslices(sec);
            z=1-linspace(0+0.5*dz, 1-0.5*dz, nslices(sec));
            
            % Normalized x positions in periodic region
            x=linspace(-0.5, 0.5, Nx);
            
            for i = 1:nslices(sec)
                
                % Compare relief with current z value
                trans = sign(relief-z(i));
                jloc = findchangepts(trans,'MaxNumChanges',Nx);
                jloc = [jloc-1;jloc];
                jno = length(jloc);
                if jno > 0
                    xtemp = x(jloc)';
                    zetas = (xtemp*[1;1]/2).';
                else
                    zetas = [-0.5,0.5];
                end
                
                % Pad edges wth \pm 0.5
                if(zetas(1) ~= -0.5)
                    zetas = [-0.5, zetas];
                end
                if(zetas(end) ~= 0.5)
                    zetas = [zetas, 0.5];
                end
                zetalist{i+loc} = zetas';
                maxlength = max(maxlength, length(zetas));
            end
        else
            error('Grating method not defined');
        end
        periodiccount = periodiccount + 1;
    else
        for i = 1:nslices(sec)
            zetalist{i+loc} = [-0.5, 0.5]';
           
        end
    end
    loc = loc + nslices(sec);
    
end

% Pad so all slices are the same size by adding 0.5 to the end of
% everything that is too short
for i = 1:Nz
    current = zetalist{i};
    len = length(current);
    if len < maxlength
        zetalist{i} = padarray(current, maxlength - len, 0.5,  'post');
    end
end

% Change cell to matrix and save to sim
sim.z.zetalist = cell2mat(zetalist);
A=1;
kappa=2;
alpha=1;
phi=0.75;
 Eg0=1.51;
 C=A*(1.62-Eg0);
VEg0 = [0,0,0,0,0,1,0,0,0];
VEg1 =strcat('0&0&0&0&0&@(x)',num2str(C),'*(0.5*(1+sin(',num2str(kappa),'*2*pi*x -',num2str(phi),'*2*pi))).^',num2str(alpha),'&0&0&0');
disp(VEg1);
sim.z.Eg = 0*sim.setup.nmz;

for i = 1:nsec
    
    tempEg0 = VEg0;
    
    if ischar(tempEg0)
        tempEg0 = strtrim(strsplit(tempEg0,'&'));
    end
    %   if sim.material.elecsec(i)==1
    Eg0 = tempEg0(i);
     
    tempEg1 = VEg1;
    if ischar(tempEg1)
        tempEg1 = strtrim(strsplit(tempEg1,'&'));
    end
    %   if sim.material.elecsec(i)==1
    Eg1 = tempEg1(i);
    
    if iscell(Eg1)
        Eg1 = Eg1{1};
        
        if Eg1(1) == '@'
            Eg1 = str2func(Eg1);
            nmz = sim.setup.nmz(matcat==i);
            if i>1
                z0 = sum(nmSec(1:i-1));
            else
                z0=0;
            end
            z = (nmz - z0)/(nmSec(i));
            
            sim.z.Eg(matcat == i) =  Eg1(z) + Eg0;
        else
            Eg1 = str2double(Eg1);
            sim.z.Eg(matcat == i) = 0*sim.z.Eg(matcat == i) +  Eg0 + Eg1;
        end
    else
        sim.z.Eg(matcat == i) = 0*sim.z.Eg(matcat == i) +  Eg0 + Eg1;
        %     end
    end
end


sim.zx.Eg = sim.z.Eg * ones(1,Nx);

nmlambda0=301;

nmlambda1=1200;
nlambda= 1 + floor((1200-301)/5);

material='Air & MgF2 & AZO &ZnO& CdS & CIGS & AZO &AZO& Silver';
% Temporary variables for inputted primary material data
material = strsplit(material,'&');
material2='Silver';
eps1 = zeros(length(sim.setup.nmz), nlambda);

% Temporary variables for inputted secondary material data
eps2 = eps1;
material2 = strsplit(material2,'&');
periodicregion = 1;
nmlambda = linspace(nmlambda0, nmlambda1, nlambda);
if nlambda ~= 1
    nmdlambda = (nmlambda1-nmlambda0)/(nlambda-1);
else
    nmdlambda = 1;
end


for i = 1:nsec
    
    % Get the material tag
    secmat = strtrim(material{i});
   
    % If flag is zero, use manually inputted data
    if strcmpi(secmat,'User')
        eps1(matcat == i, :) = eps(i);
         
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
    end
    
     
        % Check permitivity field exists
        if isfield(loadedMat,'eps')
            eps = loadedMat.eps;
        if isa(eps, 'double') && length(eps) == 1 
                % If a single permittivity then use this for all wavelenghts
                eps1(sim.setup.matcat == i, :) = eps;
            elseif isa(eps, 'double')
                % If an array of permittivities at a range of wavelengths,
                % interpolate for required wavelengths
                eps = interp1(eps(:,1), eps(:,2), nmlambda);
                eps = ones(sum(matcat == i), 1) * eps;
                eps1(matcat == i, :) = eps;
                
  elseif isa(eps, 'char')
                % If function, assume of form f(nmlambda, VEg)
                eps = str2func(eps);
                 %disp(sim.z.Eg(matcat==i));
                try 
                    eps = eps(nmlambda, sim.z.Eg(matcat==i));
                 
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























AA=0.2;
kappa=1;
alpha=5;
phi=0.75;
chi11=4.5;
sim.material = struct(...
    'nsec', 9, ...         % Number of independent section in study
    'elecsec', [0,0,0,0,1,1,0,0,0], ... % Flags for which sections to include in the electrical simulation
    'periodicsec', [0,0,0,0,0,0,0,1,0], ... % ??? grating? Flags for which sections to include in the electrical simulation
    'nmSec', [100,110,100, 70, 50, 100,10, 150, 100], ... % Lengths of the sections
    'icm3ND', [0,0,0,0,5e17,-2e16,0,0,0], ... % Doping of the sections
    ... % Material properties will be overwritten if 'material' not set to 0. If material has a range of possible bandgaps, will tailor it to specified Eg.
    'material', 'Air & MgF2 & AZO &ZnO& CdS & CIGS & AZO &AZO& Silver' , ... % Loads a material from the database, 0 for user input
    'VEg0', [0,0,0,0,0,0.941,0,0,0], ... % Baseline bandgap
    'VEg1', strcat('0&0&0&0&0&@(x)',num2str(AA),'*(0.5*(1+sin(',num2str(kappa),'*2*pi*x - ',num2str(phi),'*2*pi))).^',num2str(alpha),'&0&0&0'),... 
    ... % Electrical Properties
    'icm3Nc', 2.5e20, ... % Conduction band DOS
    'icm3Nv', 2.5e20, ... % Valence band DOS
    'cm2iVismun', 15, ...  % Electron mobility
    'cm2iVismup',  15, ...  % Hole mobility
    'VChi',4.5,...  % Electron Affinity
    'epsdc', 11.68, ... % DC relative permittivity of bulk silicon
    ... % Recombination Parameters
    'icm3isalpha', 0e-12, ... % Radiative recombination constant 1/(cm^3 s), related to beta via beta = alpha/ni^2
    'istausrhn', 2e-9, ... % SRH recombination constant 1/s
    'istausrhp', 1e-6, ... % SRH recombination constant 1/s
    'Cn', 0e-20, ... % Auger recombination constant
    'Cp', 0e-20, ... % Auger recombination constant
    'ET', 0.5, ... % Relative position of trap level in bandgap
    ... % Optical Properties
    'eps', [1, 1,1,1,1,15+2i,1, 1, 1] ... %[1+1e-9i, 4, 11.559+0.204i, 11.559+0.204i, -10.4356+0.8294i] ... % Optical Permittivity
    );

sim.optical.periodic = struct(...
    ... % Takes a vector of inputs. One entry per periodic region specified above
    'zr_flag', [1], ...                % Flag to toggle between 1 zeta, and 2 relief
    'zeta', 0.5, ...                   % Graing Duty Cycle - input from 0 (bottom) to 1 (top), output 0 to 1.
    'relief', '@(zeta) cos(2*pi*zeta)', ... % Function describing the surface relief, input -1/2 to 1/2 (I think), output trimmed to [0,1]
    'material2', 'Silver', ...         % Grating material ( 1 for silver, const otherwise)
    'eps2', -10.4356+0.8294i ...       % Grating metal permittivity
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
sim.phys.Enormconst = sqrt(1/(sim.phys.CiVimeps0 * sim.phys.misc)); % sqrt(2*sim.phys.Oeta0); % Scale E by this and H by 1/this such that incident power density is 1 W/m^2

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

% Other interal scaling factors
sim.phys.cmLs = cm_from_nm(sim.phys.nmLs); % Ls in cm
sim.phys.cm2isDs = sim.phys.cm2iVismus * sim.phys.VVs; %in cm^2s^-1
sim.phys.mAicm2Js = 1000 * sim.phys.Cq * sim.phys.VVs * sim.phys.cm2iVismus * sim.phys.icm3Ns/sim.phys.cmLs; % in mA cm^-2
sim.phys.icm2isGRs = sim.phys.VVs * sim.phys.cm2iVismus * sim.phys.icm3Ns/sim.phys.cmLs^2; % in cm^-3s^-1
sim.phys.lambda2s = sim.phys.CiVicmeps0 * sim.phys.VVth/(sim.phys.cmLs^2 * sim.phys.Cq * sim.phys.icm3Ns); % Dimensionless (needs multiplying by eps)
sim.phys.Qscale = sim.phys.Oeta0 * sim.phys.CiVimeps0/(sim.phys.Jshbar * sim.phys.Enormconst^2);









% Design the simulation

sim.optical.setup = struct(...
     ... % Optical Parameters
    'optical_toggle', 1, ... % Toggle: 1 to run optical simulation, 0 to skip
    'parforArg', inf, ...   % Numer of parallel agents (ignore)
    'epsmethod', 2, ...     % Formulaton used for eps in RCWA (2 For optimal, 1 for F(eps) and 0 (r) for invF(1/eps))
    'iepsmethod', 2, ...    % Formulaton used for ieps in RCWA (2 For optimal, 1 (r) for F(ieps) and 0 for invF(eps))
    'Ez_A_or_iE', 0, ...    % Choice for rebuilding Ez from Hy, 0 for A, 1 for iE
    'pol', 0, ...               % 0: Unpolarized, 1: p-polarized, 2: s-polarized
    ... % Analysis Options
    'RebuildFields', 0, ... % Use ifft to rebuild fields (overridden if Eplot or Hplot toggled)
    ... % Display options
    'optplot', 0, ...         % Produce plots (0: no, 1: yes)
    'epsplot', 0, ... Plot permittivity
    'Eplot', 0, ... Plot Efield
    'Hplot', 0, ... Plot Hfield
    ... % Mesh for angle and wavelength
    'degtheta0', 0, ...     % Minimum angle of Incident Light (not implemented)
    'degtheta1', 0, ...     % Maximum angle of Incident Light (not implemented)
    'ntheta', 1, ...        % Number of angles (not implemented)
    'nmlambda0', 301, ...   % Minimum wavelength in nm
    'nmlambda1', 1200, ...  % Maximum wavelength in nm
    'nlambda', 1 + floor((1200-301)/5), ...     % Number of wavelengths
    'Nx', 100, ...       % Number sample points in x-direction
    'nslices', [1, 1,1,1, 10, 50, 1,50, 1], ... % Number of slices in each sec (DesignJunction.m), uses 1 if not set.
    ... % Parameters for Nt convergence
    'CheckNtConvergence', 0, ... Use adaptive Nt (Use at own risk)
    'Nttol', 1e9,...      % Required tolerance in mA/(nm cm^2) for increasing Nt by 2
    'Nttol_lower', 1e-6, ...%  Required tolerance in mA/(nm cm^2) for decreasing Nt by 2 (used in serial calcs)
    'minNt', 8, ...         % Start value for Nt (tom:40)
    'maxNt', 12, ...        % Largest value for Nt (tom:40)
    'nosuccreq', 2, ... % Number of successes required
    'faillimit', 10, ... % Number of times worse before fail
    'ptol', 0.1 ... % Predictive stepping toler
);
    
sim.electrical.setup = struct(...
    ... % Electrical Parameters
    'electrical_toggle', 1, ... % Toggle: 1 to run electrical simulation, 0 to skip
    'pdeg', 6, ... % Polynomial degree
    'VdV', 0.01, ... % Change in volrage 
    'rerunbackward', 0, ... % Recalculate JV backward to check convergence via hystersis
    ... % Simulation meshing
    'nmdz', 5, ... % Approx. resolution of mesh
    'nx_sec', [10,50], ... % Number of mesh points in each section
    'meshfn' , '1', ...& sin(pi*x)+0.01 , ... % Function to specify relative lengths of each segment, scales to match Lz. 1 for uniform
    ... % Global polys. The poly ax^2+bx+c -> [a, b, c], for example.
    'upwind', 1, ...  % If toggled, will ignore value of t and upwind
    'tensor', 0, ... % Uses tensors internally
    't_1', 1e-6, ...
    't_2', 1e3, ...
    'quadtoggle', 1, ... Toggle quad vs lobatto integration
    'intdeg', 10, ... % Integration accuracy
    ... % Result plots
    'SolverSteps', 0, ... % Display solver steps
    'plotsolve', 0, ... % Number of steps to plot results. 0 for no output. Inf for first and last.
    'ExportSCdetails', 0, ... Export the details (bandgaps, densities, recombination etc.) at short-circuit (untested)
    'nx_plot', 2, ...
    'plot_densities', 1, ...
    'plot_fields', 1, ...
    'plot_levels', 1, ...
    'plot_currents', 1, ...
    'plot_U', 1, ...
    'plot_tau', 1, ...
  ... % Tolerance
    'dampmin', 1, ... % Set the initial damping level
    'damprelax', 1.1,... % Rate of relaxation of damping
    'trydamping', 40, ... % >=1 Number of times damping is increases to aid convergence
    'dampingincrease', 5, ... Amount to increase damping by each attempt
    'mesh_increase', 1.5, ... mesh multiplier on fail
    'refine_on_negative', 1, ...
    'refine_on_Jnoise', 1, ...
    'maxJnoiserefinefrac', 0.2, ...
    'refine_always', 0, ...
    'maxloop', 40, ... % Maximum number of iteration for convergence
    'maxpdeg', 8, ... %
    'maxtime', 60*60, ... % Maximum time allowed for the initial HDG convergence
    'mindV', 1e-6, ...
    'reltol', 1e-2 ,... % Required relative tolerance
    'abstol', 1e-4, ... % Required absolute tolerance   
    'warn', 1, ... % Display warning messages
    'Jtol', 0.1, ... Standard is 0.1. Larger values allow more noisy J profile.
    'forcepositive', 0, ... % Force the carriers to be positive
    'test_J', 0, ... %Test the Jacobian
    'PVrefinemax', 100, ... % Max refinement steps allowed
    'PVtol', 1e-3 ... % Tolerance for finding Pmax
);


%%%%%%%%%%%%%%%%%%%%%%
% Do not edit below this box %
%%%%%%%%%%%%%%%%%%%%%%

% Create timer
sim.electrical.setup.start_time = tic;
sim.electrical.setup.total_time = toc(sim.electrical.setup.start_time);

% Make sure nslices is correct lenght
diff =  sim.material.nsec - length(sim.optical.setup.nslices);
if diff > 0
    sim.optical.setup.nslices = padarray(sim.optical.setup.nslices, [0,1], diff, 'post');
end

% Set Nt to minNt
sim.optical.setup.Nt = sim.optical.setup.minNt;

% Total step counter
sim.electrical.setup.counter = 0;



sim.material.VEg0 = 1.2;
        xxeg0=sum(sim.material.VEg0);
       
        xeg0=abs((xxeg0-0.94)/(0.68));
        xeg1=abs((xxeg0-AA-0.94)/(0.68));
        xchi0=4.5-0.6*xeg0;
        xchi1=4.5-0.6*xeg1;
        xdeltachi=xchi1-xchi0;
       
       VChi0=[0,0,0,0,0,0,0,0,0];
       sim.material.VChi0=VChi0;
        sim.material.VChi1=strcat('0&0&0&0&0&@(x)',num2str(xchi0),'-',num2str(xdeltachi),'*x','&0&0&0');
      
        disp(sim.material.VChi1)

sim.material.VChi0=[0,0,0,0,0,0,0,0];

nmdz=5;
nmSec=[100,110,100, 70, 50, 100,10, 150, 100];
nslices =[1, 1,1,1, 10, 50, 1,50, 1]; 
nslices=nslices.*abs(sign(nmSec));
Nz=sum(nslices);
nmz=zeros(Nz,1);
sim.setup.nmz = nmz;
sim.setup.nmdz = nmdz;
sim.setup.Nz = Nz;
matcat=zeros(Nz,1);
sim.setup.matcat = matcat;




chi=[1,2,3,nan,nan];
disp(chi);
chi(isnan(chi))=[170.1997];
disp(chi)




