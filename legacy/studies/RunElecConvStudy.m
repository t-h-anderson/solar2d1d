% Solve the DD equations for a PN junction
% Types of variable
%   - Dimensioned (prepend with units)
%       - Defined on a section (length sec)
%       - Defined globally (length 1)
%   - Dimensionless
%       - Defined on a section (length sec)
%       - Defined globally (length 1)
%       - Defined as constant on each section (_c)
%       - Internal poly coefficient list
%       - Internal poly np matrix (_np)
%       - Standard poly np matrix (_snp)
% PN contains the junction to be simulated - PN in constan t
% sim contains the variables and simulation setup
clear();

normdeg = 1;

addpath(genpath('./'));

% Creat the default simulations


sim = DesignJunctionElecTest();
sim = DesignSimElecTest(sim);



%newsim= sim;
sim = BuildElecInput(sim);
results = 0;
for type = 8:8
    
    sim = DesignJunctionElecTest();
    sim = DesignSimElecTest(sim);

    %newsim= sim;
    sim = BuildElecInput(sim);
    
    sim.setup.start_optical = tic;
    sim.setup.start_electrical = tic;
    
    noloops = 20;
    param = [];
    clear results;
    
    switch type
        case 1
            noloops = 40;
        case 2
            noloops = 0;
        case 3
            noloops = 0;
        case 4
            noloops = 10;
        case 5
            noloops = 7;
        case 6
            noloops = 1;
        case 7
            noloops = 8;
        case 8
            noloops = 12;
        case 9
            noloops = 10;
    end
    
    
    for loopno = 1:noloops
        loopno
        % Resolution
        switch type
            case 1
                h = 20/loopno;
                sim.electrical.setup.nx_sec = ceil(sim.electrical.material.nmSec/h); % Order 0.5
                h = mean(sim.electrical.material.nmSec./sim.electrical.setup.nx_sec);
                param(loopno) = h;
            case 2
               
            case 3
               
            case 4
                sim.electrical.setup.t_1 = 10^(1-loopno); % order -4 to -5
                param(loopno) =  log10(sim.electrical.setup.t_1);
            case 5
                sim.electrical.setup.t_1 = 10^(loopno-1); % order 1.6
                param(loopno) =  log10(sim.electrical.setup.t_1);
            case 6
                sim.electrical.setup.t_1= 10^(loopno-noloops/2); % Shows convergence it to different (but similar) values
                param(loopno) =  log10(sim.electrical.setup.t_1);
            case 7
                sim.electrical.setup.pdeg = 1+loopno;
                param(loopno) = sim.electrical.setup.pdeg; % Order 2 (saturates quickly)
            case 8
                sim.electrical.setup.quadtoggle = 1;% =1 use quadrature
                sim.electrical.setup.intdeg = loopno;
                param(loopno) = sim.electrical.setup.intdeg; % Order
            case 9
                sim.electrical.setup.VdV = exp(-loopno);
                param(loopno) = -log(sim.electrical.setup.VdV); % Order 0.7 to 3 (probably depends on no. steps to refine max)
        end
        
        sim.electrical.setup.start_electrical = tic; % Start the timer. This was missing from the deprecated version.
        
        sim.electrical = RunElectrical(sim.electrical);
        
        results(loopno) = sim.electrical.results; % saave the results
        results2{loopno} = sim.electrical.coeffs; % save the coefficient vector
        %
        %sim.electrical.setup.nx_sec = [10,20,10]; % Order 0.5
        %sim.electrical.setup.large_t = 10^(1-loopno); % order -4 to -5
        %sim.electrical.setup.pdeg = 3;
        %sim.electrical.setup.intdeg = 5; % Order 3
        %sim.electrical.setup.quadtoggle = 0;
        %sim.electrical.setup.VdV = 0.01;
    end
    
    
    best = size(results,2);
    Jscerror = abs([results.mAicm2Jsc]-results(best).mAicm2Jsc)/results(best).mAicm2Jsc;
    FFerror = abs([results.FF]-results(best).FF)/results(best).FF;
    Pmaxerror = abs([results.Wim2Pmax]-results(best).Wim2Pmax)/results(best).Wim2Pmax;
    VOCerror = abs([results.VVOC]-results(best).VVOC)/results(best).VVOC;
    
    diff = length(param) - length(Jscerror);
    if diff~=0    
        pause(1)
    end
    
    try
    switch type
        case 1
            csvwrite('Jsc_nx1.csv', [param;Jscerror]);
            csvwrite('FF_nx1.csv', [param;FFerror]);
            csvwrite('Pmax_nx1.csv', [param;Pmaxerror]);
            csvwrite('VOC_nx1.csv', [param;VOCerror]);
        case 2
            csvwrite('Jsc_nx2.csv', [param;Jscerror]);
            csvwrite('FF_nx2.csv', [param;FFerror]);
            csvwrite('Pmax_nx2.csv', [param;Pmaxerror]);
            csvwrite('VOC_nx2.csv', [param;VOCerror]);
        case 3
            csvwrite('Jsc_nx3.csv', [param;Jscerror]);
            csvwrite('FF_nx3.csv', [param;FFerror]);
            csvwrite('Pmax_nx3.csv', [param;Pmaxerror]);
            csvwrite('VOC_nx3.csv', [param;VOCerror]);
        case 4
            csvwrite('Jsc_t1.csv', [param;Jscerror]);
            csvwrite('FF_t1.csv', [param;FFerror]);
            csvwrite('Pmax_t1.csv', [param;Pmaxerror]);
            csvwrite('VOC_t1.csv', [param;VOCerror]);
        case 5
            csvwrite('Jsc_t2.csv', [param;Jscerror]);
            csvwrite('FF_t2.csv', [param;FFerror]);
            csvwrite('Pmax_t2.csv', [param;Pmaxerror]);
            csvwrite('VOC_t2.csv', [param;VOCerror]);
        case 6
            csvwrite('Jsc_t3.csv', [param;Jscerror]);
            csvwrite('FF_t3.csv', [param;FFerror]);
            csvwrite('Pmax_t3.csv', [param;Pmaxerror]);
            csvwrite('VOC_t3.csv', [param;VOCerror]);
        case 7
            csvwrite('Jsc_pdeg.csv', [param;Jscerror]);
            csvwrite('FF_pdeg.csv', [param;FFerror]);
            csvwrite('Pmax_pdeg.csv', [param;Pmaxerror]);
            csvwrite('VOC_pdeg.csv', [param;VOCerror]);
        case 8
            csvwrite('Jsc_int.csv', [param;Jscerror]);
            csvwrite('FF_int.csv', [param;FFerror]);
            csvwrite('Pmax_int.csv', [param;Pmaxerror]);
            csvwrite('VOC_int.csv', [param;VOCerror]);
        case 9
            csvwrite('Jsc_dV.csv', [param;Jscerror]);
            csvwrite('FF_dV.csv', [param;FFerror]);
            csvwrite('Pmax_dV.csv', [param;Pmaxerror]);
            csvwrite('VOC_dV.csv', [param;VOCerror]);
    end
    catch
        pause();
    end
    
end

if 1
    figure()
    semilogy(param, Pmaxerror,'.');
    hold on
    semilogy(param, Jscerror,'.');
    semilogy(param, VOCerror,'.');
    semilogy(param, FFerror,'.');
    xlabel('Int_{deg}')
    ylabel('|\xi_{Int_{deg}} - \xi_{Int_{max}}|/\xi_{Int_{max}}')
    
    
    chopl = 0;
    chopr = 10;
    XP = [ones(size(Pmaxerror,2)-(chopl+chopr), 1), (param(chopl+1:end-chopr)')];
    YP = log(Pmaxerror(chopl+1:end-chopr)');
    CP = XP\YP;
    reg_P = exp(CP(1) + param * CP(2));
    semilogy(param, reg_P,'Color',[66, 134, 244]/256);
    
    XJ = [ones(size(Jscerror,2)-(chopl+chopr), 1), (param(chopl+1:end-chopr)')];
    YJ = log(Jscerror(chopl+1:end-chopr)');
    CJ = XJ\YJ;
    reg_J = exp(CJ(1) + param * CJ(2));
    semilogy(param, reg_J,'Color',[229, 18, 46]/256);
    
    XV = [ones(size(VOCerror,2)-(chopl+chopr), 1), (param(chopl+1:end-chopr)')];
    YV = log(VOCerror(chopl+1:end-chopr)');
    CV = XV\YV;
    reg_V = exp(CV(1) + param * CV(2));
    semilogy(param, reg_V,'Color',[229, 172, 18]/256);
    
    XF = [ones(size(FFerror,2)-(chopl+chopr), 1), (param(chopl+1:end-chopr)')];
    YF = log(FFerror(chopl+1:end-chopr)');
    CF = XF\YF;
    reg_F = exp(CF(1) + param * CF(2));
    semilogy(param, reg_F,'Color',[216, 28, 219]/256);
    
    
    legend('Pmax','Jsc','VOC','FF',...
        strcat('e = O(','e^{',num2str(CP(2),2),'p})'),...
        strcat('e = O(','e^{',num2str(CJ(2),2),'p})'),...
        strcat('e = O(','e^{',num2str(CV(2),2),'p})'),...
        strcat('e = O(','e^{',num2str(CF(2),2),'p})'),...
        'location','bestoutside')
else
    
    figure()
    loglog(param, Pmaxerror,'.');
    hold on
    loglog(param, Jscerror,'.');
    loglog(param, VOCerror,'.');
    loglog(param, FFerror,'.');
    
    xlabel('N_{p,x}')
    ylabel('|\xi_{N_{p,x}} - \xi_{N_{max}}|/\xi_{N_{max}}')
    
    
    chopl = 10;
    chopr = 15;
    XP = [ones(size(Pmaxerror,2)-(chopl+chopr), 1), log(param(chopl+1:end-chopr)')];
    YP = log(Pmaxerror(chopl+1:end-chopr)');
    CP = XP\YP;
    reg_P = exp(CP(1)) * param .^ CP(2);
    loglog(param, reg_P,'Color',[66, 134, 244]/256);
    
    XJ = [ones(size(Jscerror,2)-(chopl+chopr), 1), log(param(chopl+1:end-chopr)')];
    YJ = log(Jscerror(chopl+1:end-chopr)');
    CJ = XJ\YJ;
    reg_J = exp(CJ(1)) * param .^ CJ(2);
    loglog(param, reg_J,'Color',[0.8500, 0.3250, 0.0980]);
    
    XV = [ones(size(VOCerror,2)-(chopl+chopr), 1), log(param(chopl+1:end-chopr)')];
    YV = log(VOCerror(chopl+1:end-chopr)');
    CV = XV\YV;
    reg_V = exp(CV(1)) * param .^ CV(2);
    loglog(param, reg_V,'Color',[0.8500, 0.3250, 0.0980]);
    
    XF = [ones(size(FFerror,2)-(chopl+chopr), 1), log(param(chopl+1:end-chopr)')];
    YF = log(FFerror(chopl+1:end-chopr)');
    CF = XF\YF;
    reg_F = exp(CF(1)) * param .^ CF(2);
    loglog(param, reg_F,'Color',[0.8500, 0.3250, 0.0980]);
    
    legend('Pmax','Jsc','VOC','FF',...
        strcat('e = O(','N_{i,x}^{',num2str(CP(2),2),'p})'),...
        strcat('e = O(','N_{i,x}^{',num2str(CJ(2),2),'p})'),...
        strcat('e = O(','N_{i,x}^{',num2str(CV(2),2),'p})'),...
        strcat('e = O(','N_{i,x}^{',num2str(CF(2),2),'p})'),...
        'location','bestoutside')
end