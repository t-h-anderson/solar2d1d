function [ sim, cell ] = set_t(sim, cell)

t_1 = sim.setup.t_1;
t_2 = sim.setup.t_2;

% Attempt at upwinding
if sim.setup.upwind
    % error('set_t: upwinding needs to include build in field')
    
    E = coeffs2sbp(sim.coeffs.E, sim);
    n = coeffs2sbp(sim.coeffs.n,sim);
    p = coeffs2sbp(sim.coeffs.p,sim);
    lambda2 = sim.sbp.lambda2;
    %  En = sim.sbp.En;
    %  Ep = sim.sbp.Ep;
    
    taunL = zeros(sim.setup.nx, 1);
    taunL(1) = t_1;
    
    taunR = zeros(sim.setup.nx, 1);
    taunR(end) = t_1;
    
    taupL = zeros(sim.setup.nx, 1);
    taupL(1) = t_1;
    
    taupR = zeros(sim.setup.nx, 1);
    taupR(end) = t_1;
    
    PL = zeros(sim.setup.nx, 1);
    PR = zeros(sim.setup.nx, 1);
    
    EL = zeros(sim.setup.nx, 1);
    ER = zeros(sim.setup.nx, 1);
    
    for i = 1:sim.setup.nx
        
        EL(i) = polyval(E(:,i), -1);
        ER(i) = polyval(E(:,i), 1);
        
        nL0 = polyval(n(:,i),-1);
        nL1 = polyval(n(:,i),-0.9999);
        dndxL = (nL1-nL0)/0.0001;
        lambda2L = polyval(lambda2(:,i),-1);
        if dndxL ~=0
            PL(i) = abs(nL0 * EL(i)/(lambda2L * dndxL));
        end
        
        pR0 = polyval(p(:,i),0.9999);
        pR1 = polyval(p(:,i),1);
        dndxR = (pR1-pR0)/0.0001;
        lambda2R = polyval(lambda2(:,i),1);
        if dndxR ~=0
            PR(i) = abs(pR1 * ER(i)/(lambda2R * dndxR));
        end
        
    end
    
    
    taun0 = t_2;
    taun1 = t_1;
    
    taup0 = t_2;
    taup1 = t_1;
    
    for i = 1:sim.setup.nx-1
        if EL(i+1) > 0 && ER(i) > 0 && (PR(i) > 2 || PL(i+1)>2)
            taunR(i) = taun1;
            taunL(i+1) = taun0;
            taupR(i) = taup0;
            taupL(i+1) = taup1;
        elseif EL(i+1) < 0  && ER(i) < 0  && (PR(i) > 2 || PL(i+1)>2)
            taunR(i) = taun0;
            taunL(i+1) = taun1;
            taupR(i) = taup1;
            taupL(i+1) = taup0;
        else
            taunR(i) = 1;
            taunL(i+1) = 1;
            taupR(i) = 1;
            taupL(i+1) = 1;
        end
    end
    
    taunL(1) = taunL(2);
    taupL(1) = taupL(2);
    taunR(end) = taunR(end-1);
    taupR(end) = taupR(end-1);
    
    
    
    sim.global.taunl = taunL;
    sim.global.taunr = taunR;
    
    sim.global.taupl = taupL;
    sim.global.taupr = taupR;
    
    sim.global.tauphil = ones(sim.setup.nx, 1);
    sim.global.tauphir = ones(sim.setup.nx, 1);
    
else
    if ~isfield(sim.global,'taunl') || ~isfield(sim.global,'taunr') ...
            ||  ~isfield(sim.global,'taupl')  ||  ~isfield(sim.global,'taupr') ...
            ||  ~isfield(sim.global,'tauphil')  || ~isfield(sim.global,'tauphir')
        sim.global.taunl = ones(sim.setup.nx, 1);
        sim.global.taunr = ones(sim.setup.nx, 1);
        sim.global.taupl = ones(sim.setup.nx, 1);
        sim.global.taupr = ones(sim.setup.nx, 1);
        sim.global.tauphil = ones(sim.setup.nx, 1);
        sim.global.tauphir = ones(sim.setup.nx, 1);
    end
end


    sim.global.taun0 = sim.global.taunl(1);
    sim.global.taun1 = sim.global.taunr(end);
     
    sim.global.taup0 = sim.global.taupl(1);
    sim.global.taup1 = sim.global.taupr(end);
    
    sim.global.tauphi0 = sim.global.tauphil(1);
    sim.global.tauphi1 = sim.global.tauphir(end);
    
end

