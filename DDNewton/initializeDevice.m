
function [ sim ] = initializeDevice( sim )
% sim = initializeDevice(sim)
% -
% Create dimensionless version of the input and map to piecewise
% polynomials in the standard basis or piecewise constants _c. Also create other device properties.
% -
% Convert the inputs into piecewise polynomials/constants
%disp('InitialiseDevice: Only take section wise inputs, not functions so far');
% -
% Produce the matrixes which depend on the input


% If damping not set, initialise to dampmin
if ~isfield(sim.setup,'damping')
    sim.setup.damping = sim.setup.dampmin;
end

if ~isfield(sim, 'results')
    sim.results.loopcount = 1;
end
% Convert single inputs to a constant input for each section
sim.material = convert12nsec(sim.material, sim.material.nsec);

% Create polynomial basis for simulation
sim.setup.poly = lobpoly(sim.setup.pdeg);
%sim.electrical.setup.poly=rand(sim.electrical.setup.pdeg1,sim.electrical.setup.pdeg1);
%sim.electrical.setup.poly = eye(sim.electrical.setup.pdeg1);

sim.setup.ipoly = inv(sim.setup.poly);


% Simply access length of poly
sim.setup.nx = sum(sim.setup.nx_sec);
sim.setup.pdeg1 = sim.setup.pdeg+1;
sim.setup.np = sim.setup.pdeg1*sim.setup.nx;

% Calculate some scaled versions of the input
sim.setup.nmLz = sum(sim.material.nmSec);
sim.setup.Lz = sim.setup.nmLz/sim.phys.nmLs;
sim.setup.sec = sim.material.nmSec/sim.phys.nmLs;

% Create the mesh from input
% Setup mesh from input. If mesh constant, then set to uniform
% mesh on (0,Lz)
if ischar(sim.setup.meshfn)
    % If meshfun is a string then it should be a string of functions
    % separated by &, split and make into cell of functions
    meshfn_list = sim.setup.meshfn;
    meshfn_list = strtrim(strsplit(meshfn_list,'&'));
    for i=1:length(meshfn_list)
        check = meshfn_list{i};
        % If function is missing head, add it
        if ~strcmp(check(1),'@')
            meshfn_list{i} = strcat('@(x) 0*x+',check);
        end
    end
    
    % If only one function is provided, duplicate it
    if length(meshfn_list) == 1
        meshfn_backup = meshfn_list{1};
        for i = 1:length(sim.setup.nx_sec)
            meshfn_list{i} = meshfn_backup;
        end
    elseif length(meshfn_list) ~= length(sim.setup.nx)
        error('Length of meshfn does not match number of sections');
    end
else
    for i = 1:length(sim.setup.nx_sec)
        meshfn_list{i} = '@(x) 0*x+1';
    end
end

% Build the mesh
if ~isfield(sim.setup, 'refinesecs') % If not refining
    mesh = [];
    for i = 1:length(sim.setup.nx_sec)
        x = linspace(0,1,sim.setup.nx_sec(i)+1);
        x = x(1:end-1);
        meshfn = str2func(meshfn_list{i});
        mesh = [mesh, sim.material.nmSec(i) * meshfn(x)/sim.setup.nx_sec(i)];
    end
    sim.setup.mesh = mesh'/sim.phys.nmLs;
else % refine existing mesh
    if sum(sim.setup.refinesecs)>0
        count = 0;
        secno = 1;
        for i = 1:length(sim.setup.mesh)
            
            if i > sum(sim.setup.nx_sec(1:secno))
                secno = secno+1;
            end
            
            if sim.setup.refinesecs(i) == 1
                sim.setup.nx_sec(secno) = sim.setup.nx_sec(secno)+1;
                sim.setup.mesh = [sim.setup.mesh(1:i+count) ; sim.setup.mesh(i+count:end)];
                sim.setup.mesh(i+count:i+count+1) = sim.setup.mesh(i+count:i+count+1)/2;
                count = count + 1;
            end
        end
        
        
        sim.setup.nx = sum(sim.setup.nx_sec);
        sim.setup.np = sim.setup.pdeg1*sim.setup.nx;
        sim.setup.refinesecs = 0;
    end
end

% Create x0 and x1 - left and right ends of the cells.
sim.setup.x1 = zeros(length(sim.setup.mesh),1);
total = 0;
for i=1:length(sim.setup.mesh)
    total = total + sim.setup.mesh(i);
    sim.setup.x1(i) = total;
end
sim.setup.x0 = [0; sim.setup.x1(1:end-1)];



%Produce dimensionless versions
sim.setup.Vext = sim.input.VVext/ sim.phys.VVs;

% Load the data from data file if needed
sim = LoadMaterialData(sim);

sim.sbp.ND = sim.setup.damping * input2sbp(sim.material, 'icm3ND', sim.phys.icm3Ns, sim);
sim.sbp.Nc = input2sbp(sim.material, 'icm3Nc', sim.phys.icm3Ns, sim);
sim.sbp.Nv = input2sbp(sim.material, 'icm3Nv', sim.phys.icm3Ns, sim);
sim.sbp.mun = input2sbp(sim.material, 'cm2iVismun', sim.phys.cm2iVismus, sim);
sim.sbp.mup = input2sbp(sim.material, 'cm2iVismup', sim.phys.cm2iVismus, sim);

% Recombination parameters
sim.sbp.alpha = input2sbp(sim.material, 'icm3isalpha', sim.phys.icm3isGs, sim);
sim.sbp.tausrhn0 = input2sbp(sim.material, 'istausrhn', sim.phys.isTaus, sim);
sim.sbp.tausrhn1 = input2sbp(sim.material, 'istausrhn1', sim.phys.isTaus, sim);
sim.sbp.tausrhn=sim.sbp.tausrhn0+sim.sbp.tausrhn1;

sim.sbp.tausrhp0 = input2sbp(sim.material, 'istausrhp', sim.phys.isTaus, sim);
sim.sbp.tausrhp1 = input2sbp(sim.material, 'istausrhp1', sim.phys.isTaus, sim);
sim.sbp.tausrhp=sim.sbp.tausrhp0+sim.sbp.tausrhp1;

sim.sbp.Cn = input2sbp(sim.material, 'Cn', sim.phys.cm6isCs, sim);
sim.sbp.Cp = input2sbp(sim.material, 'Cp', sim.phys.cm6isCs, sim);

sim.sbp.ET = input2sbp(sim.material, 'ET', 1, sim);

% Produce bandgap and perturbation
sim.sbp.Eg0 = input2sbp(sim.material, 'VEg0', sim.phys.VVs, sim);
sim.sbp.Eg1 = sim.setup.damping * input2sbp(sim.material, 'VEg1', sim.phys.VVs, sim);
sim.sbp.Eg = sim.sbp.Eg0 + sim.sbp.Eg1;

sim.sbp.dEg = sbpdiff(sim.sbp.Eg, 1, sim);




% Produce affinity and perturbation

sim.sbp.Chi0 = input2sbp(sim.material,'VChi', sim.phys.VVs, sim);
sim.sbp.Chi1 = input2sbp(sim.material,'VChi1', sim.phys.VVs, sim);


sim.sbp.Chi=sim.sbp.Chi0+sim.sbp.Chi1;
%sim.sbp.Chi(isnan(sim.sbp.Chi))=[170.1997];
sim.sbp.dChi = sbpdiff(sim.sbp.Chi, 1, sim);


sim.sbp.En = sim.sbp.dChi;
sim.sbp.Ep = sim.sbp.dChi + sim.sbp.dEg;

sim.sbp.G = input2sbp(sim.input, 'icm3isG', sim.phys.icm3isGs, sim);
sim.coeffs.G = sbp2coeffs(sim.sbp.G, sim);

sim.setup.dV = sim.setup.VdV / sim.phys.VVs;


% Produce the continuous quasi Fermi-level matrices
sim.global.DeltaChi = zeros(sim.setup.nx-1);
sim.global.DeltaChiEg = zeros(sim.setup.nx-1);
%Tom please look below, I feel this is one spot where possible problem can occur 
for i = 1:sim.setup.nx-1
    sim.global.DeltaChi(i,i) =  1./exp(polyval(sim.sbp.Chi(:,i), 1) - polyval(sim.sbp.Chi(:,i+1),-1));
    sim.global.DeltaChiEg(i,i) =  1./(sim.global.DeltaChi(i,i)) * exp((polyval(sim.sbp.Eg(:,i), 1) - polyval(sim.sbp.Eg(:,i+1),-1))); 
end

% Update if they exist
if isfield(sim.global, 'SL')
    if size(sim.global.SL, 1) == size(sim.global.DeltaChi,2)
        sim.global.S_psi_n = sim.global.DeltaChi * sim.global.SL - sim.global.SR;
        sim.global.S_psi_p = sim.global.DeltaChiEg * sim.global.SL - sim.global.SR;
        sim.global.S_psi = sim.global.SL - sim.global.SR;
    end
end

% Create compound sbp which are repeatedly used in the simulation
ND = sim.sbp.ND;
%sim.coeffs.ND = sbp2coeffs(ND, sim);
sim.sbp.iND = polyfunc(@(x)1./x, ND);
sim.sbp.ND2 = sbptimes(ND, ND, sim);

NcNv = sbptimes(sim.sbp.Nc, sim.sbp.Nv, sim);
coeffpart = polyfunc(@(x)sqrt(x), NcNv);
exppart = polyfunc(@(x)exp(-x), sim.sbp.Eg/2);
sim.sbp.ni = sbptimes(coeffpart, exppart, sim);
sim.sbp.ni2 = polyfunc(@(x)x.^2, sim.sbp.ni);
sim.sbp.ini2 = polyfunc(@(x)1./x, sim.sbp.ni2);
sim.sbp.lambda2 = input2sbp(sim.material, 'epsdc', 1/sim.phys.lambda2s, sim);

sim.sbp.E1pos = sbptimes(sim.sbp.ET, sim.sbp.Eg, sim);
sim.sbp.E2pos = sim.sbp.Eg - sim.sbp.E1pos;

sim.sbp.n1 = sbptimes(polyfunc(@(x)exp(-x), sim.sbp.E1pos), sim.sbp.Nc, sim);
sim.sbp.p1 = sbptimes(polyfunc(@(x)exp(-x), sim.sbp.E2pos), sim.sbp.Nv, sim);


sim.coeffs.ni = sbp2coeffs(sim.sbp.ni, sim);
%sim.coeffs.ni2 = sbp2coeffs(sim.sbp.ni2, sim);
sim.coeffs.n1 = sbp2coeffs(sim.sbp.n1, sim);
sim.coeffs.p1 = sbp2coeffs(sim.sbp.p1, sim);

sim.coeffs.tausrhn = sbp2coeffs(sim.sbp.tausrhn, sim);
sim.coeffs.tausrhp = sbp2coeffs(sim.sbp.tausrhp, sim);

sim.coeffs.Cn = sbp2coeffs(sim.sbp.Cn, sim);
sim.coeffs.Cp = sbp2coeffs(sim.sbp.Cp, sim);

% Weighted current mass matrixes
sim.sbp.c_n = polyfunc(@(x) 1./x, sim.sbp.mun);
sim.sbp.c_p = - polyfunc(@(x) 1./x, sim.sbp.mup);
sim.sbp.c_phi = - polyfunc(@(x) 1./x, sim.sbp.lambda2);

sim.coeffs.ND = sbp2coeffs(sim.sbp.ND, sim);

end


