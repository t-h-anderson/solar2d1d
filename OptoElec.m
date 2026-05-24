
function [fom] = OptoElec(unused, params)

% Fom stores the opjective variable first, followed by a list of
% non-optimised objectives. Tip: include the objective twice so it's also
% stored in the sublist.
fom = cell(11, 1);

% Create the simulation
sim = DesignJunction();
sim = DesignSim(sim);

for i = 1:length(fom)
    fom{i} = nan;
end
optsuccessflag = 0;

% Set parameter mapping
%sim.material.nmSec(6) = params.nmLi;
disp(sim.material.nmSec(6));
sim.material.nmSec(8) = params.nmLg;
sim.optical.periodic.zeta = params.zeta;
sim.input.nmLx = params.nmLx;

% Bandgap-profile mode for the absorber section. Was a hard-coded
% `switch 3` in the original; lift it onto params so callers can pick
% explicitly without editing the source.
%   1 = linear ramp across three sub-regions
%   2 = single flat bandgap (Eg_p)
%   3 = sinusoidal profile with amplitude A, freq kappa, phase phi, exponent alpha
if isfield(params, 'bandgap_profile_mode')
    bandgap_profile_mode = params.bandgap_profile_mode;
else
    bandgap_profile_mode = 3;
end

switch bandgap_profile_mode
    case 1
        Eg0 = params.Eg0;
        Eg1 = params.Eg1;
        Egrange = [Eg0, Eg1];
        
        m = num2str((Eg1-Eg0)/(params.nmLi+40));
        f1 = strcat('@(x)', m, '*x*20');
        f2 = strcat('@(x) ',num2str(m),'*(',num2str(params.nmLi),'*x + 20)');
        f3 = strcat('@(x) ',m,'*(20*x +', num2str(params.nmLi), '+ 20)');
        
        sim.material.VEg1 = strcat('0&0&',f1,'&',f2,'&',f3,'&0&0');
        sim.material.VEg0 = Eg0;
    case 2
        Eg_p = params.Eg_p;
        Egrange = [Eg_p];
        sim.material.VEg0 = Eg_p;
       
        
    case 3
        Eg0=params.Eg0;
        A = params.A;
        AA=A*(1.49-Eg0);
        sim.material.VEg0(6) = Eg0;
        disp(sim.material.VEg0);
        xxeg0=sim.material.VEg0(6);
        
        xeg0=abs((xxeg0-0.91)/(0.58));
        kappa = params.kappa;
        phi = params.phi;
        alpha = params.alpha;
        sim.material.VEg1 = strcat('0&0&0&0&0&@(x)',num2str(AA),'*(0.5*(1+sin(',num2str(kappa),'*2*pi*x - ',num2str(phi),'*2*pi))).^',num2str(alpha),'&0&0&0');
        
        Egrange=[Eg0,Eg0+AA];
        disp(Egrange);
        
        xeg1=abs((xxeg0+AA-0.91)/(0.58));
        xchi0=4.46-0.16*xeg0;
        disp(xchi0);
        xchi1=4.46-0.16*xeg1;
        xdeltachi=xchi0-xchi1;
        sim.material.VChi1=strcat('0&0&0&0&0&@(x)',num2str(xchi0),'-',num2str(xdeltachi),'*(0.5*(1+sin(',num2str(kappa),'*2*pi*x - ',num2str(phi),'*2*pi))).^',num2str(alpha),'&0&0&0');
        
       %Bandgap dependent Defect Density%%%%%%%%%%%%%%%%%%%
       xtausrhn0=2e-9-1.9e-9*xeg0;
       xtausrhn1=2e-9-1.9e-9*xeg1;
       xdeltasrhn=xtausrhn0-xtausrhn1;
        disp(xdeltasrhn);
       sim.material.istausrhn1=strcat('0&0&0&0&0&@(x)',num2str(xtausrhn0),'-',num2str(xdeltasrhn),'*(0.5*(1+sin(',num2str(kappa),'*2*pi*x - ',num2str(phi),'*2*pi))).^',num2str(alpha),'&0&0&0');
       xtausrhp0=2e-9-1.9e-9*xeg0;
       xtausrhp1=2e-9-1.9e-9*xeg1;
       xdeltasrhp=xtausrhp0-xtausrhp1;
       disp(xtausrhp0);
       sim.material.istausrhp1=strcat('0&0&0&0&0&@(x)',num2str(xtausrhp0),'-',num2str(xdeltasrhp),'*(0.5*(1+sin(',num2str(kappa),'*2*pi*x - ',num2str(phi),'*2*pi))).^',num2str(alpha),'&0&0&0');
end

sim = BuildOptInput(sim);

% Scale mesh to keep same resolution
sim.optical.setup.nslices(8) = ceil(params.nmLg/5);
%sim.optical.setup.nslices(6) = ceil(params.nmLi/5);
%disp(sim.optical.setup.nslices(6));
%sim.optical.setup.Nx = ceil(params.nmLx/2);
%sim.electrical.setup.nx_sec(3)=ceil(params.nmLi/5);
%disp(sim.electrical.setup.nx_sec(3));
%nmdz = sim.electrical.setup.nmdz; %
%if nmdz ~= 0
 %   sim.electrical.setup.nx_sec = ceil(sim.material.nmSec/nmdz);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start optimisation
% no need to edit below here unless number of subvalues changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    sim.optical.setup.start_optical = tic;
    if sim.optical.setup.optical_toggle
        % Copy input to optical
        sim.optical = RCWAinitilisation(sim.optical);
        sim.optical = RunOptical(sim.optical);
        
        % Copy the results to the DEA object
        fom{2} = sim.optical.results.mAicm2JOpt;
        fom{3} = sim.optical.results.mAicm2JOpt_s;
        fom{4} = sim.optical.results.mAicm2JOpt_p;
        optsuccessflag = 1;
        sim.input.icm3isG = sim.optical.results.icm3isG;
        sim.input.mAicm2JOpt = sim.optical.results.mAicm2JOpt;
    end
    sim.optical.results.total_optical = toc(sim.optical.setup.start_optical);
    
    disp('Optical Done')
    
     % Run electrical model
    sim.electrical.setup.start_electrical = tic;
    sim = BuildElecInput(sim);
    if sim.electrical.setup.electrical_toggle   
        sim.electrical = RunElectrical(sim.electrical);
    end
    sim.electrical.results.total_electrical = toc(sim.electrical.setup.start_electrical);
    %sim.electrical.results
    
    disp('Electrical Done')
    
    eta = sim.electrical.results.Wim2Pmax/10;
    Jsc = sim.electrical.results.mAicm2Jsc;
    JscOpt = sim.electrical.results.JscOpt;
    Voc = sim.electrical.results.VVOC;
    FF = sim.electrical.results.FF;
    Vmax = sim.electrical.results.VVmax;
    
    fom{1} = eta;
    fom{5} = eta;
    fom{6} = Vmax;
    fom{7} = Jsc;
    fom{8} = Voc;
    fom{9} = FF;
    fom{10} = sim.optical.results.total_optical;
    fom{11} = sim.electrical.results.total_electrical;
    
    
    
    if eta < 0
        disp('eta < 0');
        error('eta < 0');
    elseif eta > 100
        disp('eta>100');
        error('eta>100');
    elseif Jsc < 0
        disp('Jsc < 0');
        error('Jsc < 0');
    elseif Jsc > JscOpt
        disp('Jsc > JscOpt');
        error('Jsc > JscOpt');
    elseif JscOpt < 0
        disp('JscOpt < 0');
        error('JscOpt < 0');
    elseif Voc < 0
        disp('Voc < 0');
        error('Voc < 0');
    elseif Voc > max(Egrange)
        disp('Voc > possible');
        error('Voc > possible');
    elseif FF <=0
        disp('FF<0');
        error('FF<0');
    elseif FF > 1
        disp('FF>1');
        error('FF>1');
    elseif isempty(FF)
        error('JV Failed')
    elseif isempty(Voc)
        error('JV Failed')
    elseif isempty(Jsc)
        error('JV Failed')
    end
    disp(sprintf('eta: %.2f%%, Voc: %.2f, Jsc:%.2f, JscOpt:%.2f, FF: %.2f \n',eta, Voc, Jsc, JscOpt, FF));
 catch me
     
     if optsuccessflag
         disp('Likely Electrical Computation Error');       
         fom{1} = nan;
         fom{5} = nan;
         fom{6} = nan;
         fom{7} = nan;
         fom{8} = nan;
         fom{9} = nan;
         fom{11} = nan;
     else
         disp('Optical Computation Error');       
         fom{1} = nan;
         fom{2} = nan;
         fom{3} = nan;
         fom{4} = nan;
         fom{5} = nan;
         fom{6} = nan;
         fom{7} = nan;
         fom{8} = nan;
         fom{9} = nan;
         fom{10} = nan;
         fom{11} = nan;
     end
     
     disp(me);
 end

end





