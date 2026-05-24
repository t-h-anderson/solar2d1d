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










