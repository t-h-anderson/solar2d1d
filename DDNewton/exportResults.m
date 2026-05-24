function exportResults(sim)
% Choose which results to display. toggle should be a list
% of binary options
% [ densities, currents, fields,
% Display running results

here = pwd;
%here = "/home/WIN/tomha/Dropbox/PV_Model/Anderson_JoCP/simdetails/";
%here = "C:\Users\Tom\Dropbox\PV_Model\Anderson_JoCP\data\simdetails\";

% Densities
    n = coeffs2sbp(sim.coeffs.n, sim);
    p = coeffs2sbp(sim.coeffs.p, sim);
    rho = coeffs2sbp(sim.coeffs.p-sim.coeffs.n+sim.coeffs.ND,sim);
    
        
    pdata = makeresults(p, sim.phys.nmLs, sim.phys.icm3Ns, sim, 10);
    csvwrite(strcat(here,'pdata.csv'), pdata);
    
    p_hatdata = makeresults(sim.coeffs.p_hat', sim.phys.nmLs, sim.phys.icm3Ns, sim, 0);
    csvwrite(strcat(here,'p_hatdata.csv'), p_hatdata);
    
    ndata = makeresults(n, sim.phys.nmLs,  sim.phys.icm3Ns, sim, 10);
    csvwrite(strcat(here,'ndata.csv'), ndata);
    
    n_hatdata = makeresults(sim.coeffs.n_hat', sim.phys.nmLs, sim.phys.icm3Ns, sim, 0);
    csvwrite(strcat(here,'n_hatdata.csv'), n_hatdata);
    
    rhodata = makeresults(rho, sim.phys.nmLs,  sim.phys.icm3Ns, sim, 10);
    csvwrite(strcat(here,'rhodata.csv'), rhodata);
    
    ni = makeresults(sim.sbp.ni, sim.phys.nmLs, sim.phys.icm3Ns, sim, 10);
    csvwrite(strcat(here,'nidata.csv'), ni);

    % Currents
    Jn = coeffs2sbp(sim.coeffs.Jn, sim);
    Jp = coeffs2sbp(sim.coeffs.Jp, sim);
        
    Jndata = makeresults(Jn, sim.phys.nmLs, sim.phys.mAicm2Js, sim, 10);    
    csvwrite(strcat(here,'Jndata.csv'), Jndata);
    
    Jpdata = makeresults(Jp, sim.phys.nmLs, sim.phys.mAicm2Js, sim, 10);
    csvwrite(strcat(here,'Jpdata.csv'), Jpdata);
  
        for i = 1:(sim.setup.nx -1)
            Jn_hatr(i) = polyval(Jn(:,i),1) + sim.global.taunr(i) * (polyval(n(:,i), 1) - sim.coeffs.n_hat(i));
            Jn_hatl(i) = polyval(Jn(:,i+1),-1) + sim.global.taunl(i) * (polyval(n(:,i+1), -1) - sim.global.DeltaChi(i,i)*sim.coeffs.n_hat(i));
            
            Jp_hatr(i) = polyval(Jp(:,i),1) + sim.global.taupr(i) * (polyval(p(:,i), 1) - sim.coeffs.p_hat(i));
            Jp_hatl(i) = polyval(Jp(:,i+1),-1) + sim.global.taupl(i) * (polyval(p(:,i+1), -1) - sim.global.DeltaChiEg(i,i)*sim.coeffs.p_hat(i));
        end
        
        Jnr_hatdata = makeresults(Jn_hatr, sim.phys.nmLs, sim.phys.mAicm2Js, sim, 0);
        csvwrite(strcat(here,'Jnr_hatdata.csv'), Jnr_hatdata);
        
        Jpr_hatdata = makeresults(Jp_hatr, sim.phys.nmLs, sim.phys.mAicm2Js, sim, 0);
        csvwrite(strcat(here,'Jpr_hatdata.csv'), Jpr_hatdata);
        
        Jnl_hatdata = makeresults(Jn_hatl, sim.phys.nmLs, sim.phys.mAicm2Js, sim, 0);
        csvwrite(strcat(here,'Jnl_hatdata.csv'), Jnl_hatdata);
        
        Jpl_hatdata = makeresults(Jp_hatl, sim.phys.nmLs, sim.phys.mAicm2Js, sim, 0);
        csvwrite(strcat(here,'Jpl_hatdata.csv'), Jpl_hatdata);
        
    % Electrics
    E = coeffs2sbp(sim.coeffs.E, sim);
    phi = coeffs2sbp(sim.coeffs.phi, sim);
    phin = phi + sim.sbp.dEg;
    phip = phin + sim.sbp.dChi;
    
    sim.setup.lp = 0;
    Edata = makeresults(E, sim.phys.nmLs, 1, sim, 10);
    csvwrite(strcat(here,'Edata.csv'), Edata);
    
    phidata = makeresults(phi, sim.phys.nmLs, sim.phys.VVs, sim, 10);
    csvwrite(strcat(here,'phidata.csv'), phidata);
    
    phindata = makeresults(phin, sim.phys.nmLs, sim.phys.VVs, sim, 10);
    csvwrite(strcat(here,'phindata.csv'), phindata);
    
    phipdata = makeresults(phip, sim.phys.nmLs, sim.phys.VVs, sim, 10);
    csvwrite(strcat(here,'phipdata.csv'), phipdata);
    
    phi_hatdata = makeresults(sim.coeffs.phi_hat', sim.phys.nmLs, sim.phys.VVs, sim, 0);
    csvwrite(strcat(here,'phi_hatdata.csv'), phi_hatdata);
    
    %plotlevels(sim, figno, 1)
    chi =  sbp2space(sim.sbp.Chi, sim)';
    x = sim.phys.nmLs * chi(1,:);
        
    Evac = -coeffs2space(sim.coeffs.phi, sim)';
    n = coeffs2space(sim.coeffs.n, sim)';
    p = coeffs2space(sim.coeffs.p, sim)';
    ni = sbp2space(sim.sbp.ni, sim)';
    Eg = sbp2space(sim.sbp.Eg, sim)';

    Nc = sbp2space(sim.sbp.Nc,sim)';
    Nv = sbp2space(sim.sbp.Nv,sim)';
    
    
    % Distance from E_fn to phi at LHS
    Vc0 = - log(sim.bc.n0/Nc(2,1));
    EF0 = Evac(2,1) - chi(2,1) - Vc0;
    
    Evac1 = [Evac(1,:); Evac(2,:)];
    Ec = [Evac(1,:); Evac(2,:) - chi(2,:)];
    Ei = [Evac(1,:); Ec(2,:) - 0.5*Eg(2,:) - 0.5*log(Nc(2,:)./Nv(2,:))];
    Fn = [Evac(1,:); Ei(2,:) + log(n(2,:)./ni(2,:))];
    Fp = [Evac(1,:); Ei(2,:) - log(p(2,:)./ni(2,:))];  
    Ev = [Evac(1,:); Ec(2,:) - Eg(2,:)];
    
    Evacdata = [x(1,:); sim.phys.VVs * Evac1(2,:)]'; % Space field  
    csvwrite(strcat(here,'Evacdata.csv'), Evacdata);

    Ecdata  = [x(1,:); sim.phys.VVs * real(Ec(2,:))]'; % Conduction Band
    csvwrite(strcat(here,'Ecdata.csv'), Ecdata); 
    
    Fndata = [x(1,:); sim.phys.VVs * real(Fn(2,:))]'; % Electron Fermi level
    csvwrite(strcat(here,'Fndata.csv'), Fndata);
    
    Eidata = [x(1,:); sim.phys.VVs * real(Ei(2,:))]'; % Intrinsic energy level
    csvwrite(strcat(here,'Eidata.csv'), Eidata);
    
    Fpdata = [x(1,:); sim.phys.VVs * real(Fp(2,:))]'; % Hole Fermi level
    csvwrite(strcat(here,'Fpdata.csv'), Fpdata);
    
    Evdata = [x(1,:); sim.phys.VVs * real(Ev(2,:))]'; % Valence Band
    csvwrite(strcat(here,'Evdata.csv'), Evdata);
    
    
    
    % Generation
    G = coeffs2sbp(sim.coeffs.G, sim);
    coeff = sbptimes(sim.sbp.ini2, sim.sbp.alpha, sim);
    n = coeffs2sbp(sim.coeffs.n, sim);
    p = coeffs2sbp(sim.coeffs.p, sim);
    np = sbptimes(n, p, sim);
    
    % Direct
    numerator = np - sim.sbp.ni2;
    RD = sim.setup.damping * sbptimes(coeff, numerator, sim);
    
    % SRH
    denominator = sbptimes(sim.sbp.tausrhp, n + sim.sbp.n1, sim) + sbptimes(sim.sbp.tausrhn, p + sim.sbp.p1, sim);
    denominator = polyfunc(@(x)1./x, denominator);
    R_SRH = sim.setup.damping * sbptimes(numerator, denominator, sim);
    
    % Auger
    CnnCpp = sbptimes(sim.sbp.Cn, n, sim) + sbptimes(sim.sbp.Cp, p, sim);
    R_Aug = sim.setup.damping * (sbptimes(CnnCpp, numerator, sim));
    
    R = RD + R_SRH + R_Aug;

    Gdata = makeresults(G, sim.phys.nmLs, sim.phys.icm3isGs, sim, 10);
    csvwrite(strcat(here,'Gdata.csv'), Gdata);
    
    Rdata = makeresults(R, sim.phys.nmLs, sim.phys.icm3isGs, sim, 10);
    csvwrite(strcat(here,'Rdata.csv'), Rdata);
    
    RDdata = makeresults(RD, sim.phys.nmLs, sim.phys.icm3isGs, sim, 10);
    csvwrite(strcat(here,'RDdata.csv'), RDdata);
    
    RSRHdata = makeresults(R_SRH, sim.phys.nmLs, sim.phys.icm3isGs, sim, 10);
    csvwrite(strcat(here,'RSRHdata.csv'), RSRHdata);
    
    RAugdata = makeresults(R_Aug, sim.phys.nmLs, sim.phys.icm3isGs, sim, 10);
    csvwrite(strcat(here,'RAugdata.csv'), RAugdata);
    
end