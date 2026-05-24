



function [sim] = BuildMaterial(sim)
% BUILDMATERIAL Categorises the material regions

%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim.setup.nslices = sim.setup.nslices .* abs(sign(sim.material.nmSec));%%%%%%%convert all thickness to unity
%sim.setup.nslices
Nz = sum(sim.setup.nslices);
nmz = zeros(Nz, 1);
nmdz = zeros(Nz, 1);
matcat = zeros(Nz, 1);

z0 = 0;
z1 = 0;
index1 = 1;

elecmask = zeros(Nz, sim.setup.Nx);
nmdzmask = zeros(Nz, sim.setup.Nx);
for sec = 1:length(sim.setup.nslices)
    
    index2 = index1 + sim.setup.nslices(sec) - 1;
    nmdz(index1:index2) = sim.material.nmSec(sec)/sim.setup.nslices(sec);%%%%%%dz in each discrete section
    
    z1 = z1 + sim.material.nmSec(sec);%%%%%%adding thicknesses
    if z1 ~= z0
        newnmz = linspace(z0 + nmdz(index1)/2, z1 - nmdz(index1)/2, sim.setup.nslices(sec));
    else
        newnmz = [];
    end
    nmz(index1:index2) = newnmz;
    matcat(index1:index2) = sec;
    
    elecmask(index1:index2, :) = sim.material.elecsec(sec);
    nmdzmask(index1:index2, :) = sim.material.nmSec(sec)/sim.setup.nslices(sec);
    
    z0=z1;
    index1 = index2+1;
end

% Create z0 coordinates
sec = 1;
sim.setup.nmz0 = 0;
while sec <= length(sim.material.elecsec) && sim.material.elecsec(sec)==0
    sim.setup.nmz0 = sim.setup.nmz0 + sim.material.nmSec(sec);
    sec = sec + 1;
end

sim.setup.nmz1 = sim.setup.nmz0;
while sec <= length(sim.material.elecsec) && sim.material.elecsec(sec)==1
    sim.setup.nmz1 = sim.setup.nmz1 + sim.material.nmSec(sec);
    sec = sec + 1;
end

sim.setup.nmLz = sim.setup.nmz1 - sim.setup.nmz0;
sim.setup.Lz = sim.setup.nmLz / sim.phys.nmLs;

sim.setup.nmz = nmz;
sim.setup.nmdz = nmdz;
sim.setup.Nz = Nz;
sim.setup.matcat = matcat;

sim.zx.elecmask = elecmask;
sim.zx.nmdzmask = nmdzmask;

sim.z.elecmask = elecmask(:,1);
sim.z.nmdzmask = nmdzmask(:,1);

zetas = linspace(-0.5,0.5,sim.setup.Nx);
zetalist = [];
periodiccount =1;
loc = 0;
maxlength = 2;

% Loop through layers
for sec = 1:sim.material.nsec
    
    % If periodic section
    if sim.material.periodicsec(sec) == 1
        
        % If def by zeta
        if sim.periodic.zr_flag(periodiccount) == 1
            
            % If number
            if isa(sim.periodic.zeta(periodiccount), 'double')
                % Set all layers to number
                for i = 1:sim.setup.nslices(sec)
                    zetas = [-0.5*sim.periodic.zeta(periodiccount), 0.5*sim.periodic.zeta(periodiccount)];
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
                zetacells = strsplit(sim.periodic.zeta,'&');
                zetaFn = zetacells{periodiccount};
                
                
                
                zetaFn = str2func(zetaFn);
                
                dz = 1/ sim.setup.nslices(sec);
                z=1-linspace(0+0.5*dz, 1-0.5*dz, sim.setup.nslices(sec));
                for i = 1:sim.setup.nslices(sec)
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
        elseif sim.periodic.zr_flag(periodiccount) == 2 % Specify sim.z.zeta via input surface
            
            % Load function for section
            reliefcells = strsplit(sim.periodic.relief,'&');
            reliefFn = reliefcells{periodiccount};
            reliefFn = str2func(reliefFn);
            zeta = linspace(-0.5, 0.5, sim.setup.Nx);
            
            % Surface g(x)
            relief = reliefFn(zeta);
            
            % Normalized z positions in periodic region
            dz = 1/ sim.setup.nslices(sec);
            z=1-linspace(0+0.5*dz, 1-0.5*dz, sim.setup.nslices(sec));
            
            % Normalized x positions in periodic region
            x=linspace(-0.5, 0.5, sim.setup.Nx);
            
            for i = 1:sim.setup.nslices(sec)
                
                % Compare relief with current z value
                trans = sign(relief-z(i));
                jloc = findchangepts(trans,'MaxNumChanges',sim.setup.Nx);
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
        for i = 1:sim.setup.nslices(sec)
            zetalist{i+loc} = [-0.5, 0.5]';
        end
    end
    loc = loc + sim.setup.nslices(sec);
    
end

% Pad so all slices are the same size by adding 0.5 to the end of
% everything that is too short
for i = 1:sim.setup.Nz
    current = zetalist{i};
    len = length(current);
    if len < maxlength
        zetalist{i} = padarray(current, maxlength - len, 0.5,  'post');
    end
end

% Change cell to matrix and save to sim
sim.z.zetalist = cell2mat(zetalist);
end

