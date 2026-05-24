



function sim = RCWAinitilisation(sim)
%Run the RCWA algorithm for a range of wavelengths to build the absorption spectrum

% Convert single inputs to a constant input for each section
sim.material = convert12nsec(sim.material, sim.material.nsec);
sim.periodic = convert12nsec(sim.periodic, sum(sim.material.periodicsec));
sim.periodic.nsec = sum(sim.material.periodicsec == 1);

% Create wavelength and wavenumber vectors
sim.setup.nmlambda = linspace(sim.setup.nmlambda0, sim.setup.nmlambda1, sim.setup.nlambda);
if sim.setup.nlambda ~= 1
    sim.setup.nmdlambda = (sim.setup.nmlambda1-sim.setup.nmlambda0)/(sim.setup.nlambda-1);
else
    sim.setup.nmdlambda = 1;
end
sim.setup.inmk0 = 2.*pi./sim.setup.nmlambda;
sim.setup.radtheta = rad_from_deg(linspace(sim.setup.degtheta0, sim.setup.degtheta1, sim.setup.ntheta));

% Create Material Vectors
sim.setup.nmx = linspace(-0.5*sim.input.nmLx, 0.5*sim.input.nmLx, sim.setup.Nx);
sim.setup.nmdx = sim.input.nmLx/(sim.setup.Nx - 1);
sim = BuildMaterial(sim);

% Calculate Eg at junction locations
sim.z.Eg = 0*sim.setup.nmz;
for i = 1:sim.material.nsec
    
    tempEg0 = sim.material.VEg0;
    if ischar(tempEg0)
        tempEg0 = strtrim(strsplit(tempEg0,'&'));
    end
    %   if sim.material.elecsec(i)==1
    Eg0 = tempEg0(i);
     
    tempEg1 = sim.material.VEg1;
    if ischar(tempEg1)
        tempEg1 = strtrim(strsplit(tempEg1,'&'));
    end
    %   if sim.material.elecsec(i)==1
    Eg1 = tempEg1(i);
    
    if iscell(Eg1)
        Eg1 = Eg1{1};
        
        if Eg1(1) == '@'
            Eg1 = str2func(Eg1);
            nmz = sim.setup.nmz(sim.setup.matcat==i);
            if i>1
                z0 = sum(sim.material.nmSec(1:i-1));
            else
                z0=0;
            end
            z = (nmz - z0)/(sim.material.nmSec(i));
            
            sim.z.Eg(sim.setup.matcat == i) =  Eg1(z) + Eg0;
        else
            Eg1 = str2double(Eg1);
            sim.z.Eg(sim.setup.matcat == i) = 0*sim.z.Eg(sim.setup.matcat == i) +  Eg0 + Eg1;
        end
    else
        sim.z.Eg(sim.setup.matcat == i) = 0*sim.z.Eg(sim.setup.matcat == i) +  Eg0 + Eg1;
        %     end
    end
end


sim.zx.Eg = sim.z.Eg * ones(1, sim.setup.Nx);








% Load permittivites and create eps(z,l). Need to be two different values
% for grating
sim = LoadEpsData(sim);

% Matrices for permittivity and fourier of permittivity
%sim.lxz.eps = zeros(sim.setup.nlambda, sim.setup.Nx, sim.setup.Nz);

sim.phys.Wim2inm2Sl = Wim2inmAM15G(sim.setup.nmlambda);

% Electric Field Variables
sim.zfl.Ex = cell(sim.setup.nlambda,1);
sim.zfl.Ey = cell(sim.setup.nlambda,1);
sim.zfl.Ez = cell(sim.setup.nlambda,1);

sim.zx.eps = zeros(sim.setup.Nz, sim.setup.Nx);%, sim.setup.nlambda);
sim.zx.Ex = zeros(sim.setup.Nz, sim.setup.Nx); %, sim.setup.nlambda);
sim.zx.Ey = zeros(sim.setup.Nz, sim.setup.Nx); %, sim.setup.nlambda);
sim.zx.Ez = zeros(sim.setup.Nz, sim.setup.Nx); %, sim.setup.nlambda);
sim.zx.Hy = zeros(sim.setup.Nz, sim.setup.Nx); %, sim.setup.nlambda);
sim.zx.Es = zeros(sim.setup.Nz, sim.setup.Nx); %, sim.setup.nlambda);
sim.zx.Ep = zeros(sim.setup.Nz, sim.setup.Nx); %, sim.setup.nlambda);

% Empty matrixes for absorption and reflection
sim.results.t = cell(sim.setup.nlambda,1);
sim.results.r = cell(sim.setup.nlambda,1);

sim.results.st = cell(sim.setup.nlambda,1);
sim.results.sr = cell(sim.setup.nlambda,1);

sim.results.pt = cell(sim.setup.nlambda,1);
sim.results.pr = cell(sim.setup.nlambda,1);

sim.results.T = zeros(sim.setup.nlambda,1);
sim.results.R = zeros(sim.setup.nlambda,1);
sim.results.A = zeros(sim.setup.nlambda, 1);
sim.results.Aelec = zeros(sim.setup.nlambda, 1);

sim.results.sT = zeros(sim.setup.nlambda,1);
sim.results.sR = zeros(sim.setup.nlambda,1);
sim.results.sA = zeros(sim.setup.nlambda, 1);
sim.results.sAelec = zeros(sim.setup.nlambda, 1);

sim.results.Ajunc = zeros(sim.setup.nlambda,1);
sim.results.Ajunc_s = zeros(sim.setup.nlambda,1);
sim.results.Ajunc_p = zeros(sim.setup.nlambda,1);
sim.results.Arest = zeros(sim.setup.nlambda,1);

sim.results.pT = zeros(sim.setup.nlambda,1);
sim.results.pR = zeros(sim.setup.nlambda,1);
sim.results.pA = zeros(sim.setup.nlambda, 1);
sim.results.sAelec = zeros(sim.setup.nlambda, 1);
end


