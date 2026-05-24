function sim = displayResults(sim)
% Choose which results to display. toggle should be a list
% of binary options
% [ densities, currents, fields,
% Display running results

figure(1);
figno = 1;

if sim.setup.plot_densities
    n = coeffs2sbp(sim.coeffs.n, sim);
    p = coeffs2sbp(sim.coeffs.p, sim);
    rho = coeffs2sbp(sim.coeffs.p-sim.coeffs.n+sim.coeffs.ND,sim);
    
    sim.setup.lp = 1;
    
    makeplot(p, 'Densities', sim.phys.nmLs, sim.phys.icm3Ns, sim, figno, 1,'blue',10);
    makeplot(sim.coeffs.p_hat', '', sim.phys.nmLs, sim.phys.icm3Ns, sim, figno,0,'b.',0);
    makeplot(n, '', sim.phys.nmLs,  sim.phys.icm3Ns, sim, figno, 0,'red',10);
    makeplot(sim.coeffs.n_hat', '', sim.phys.nmLs, sim.phys.icm3Ns, sim, figno, 0,'r.',0);
    makeplot(rho, '', sim.phys.nmLs,  sim.phys.icm3Ns, sim, figno, 0,'green',10);
    makeplot(sim.sbp.ni, '', sim.phys.nmLs, sim.phys.icm3Ns, sim, figno, 0,'black--',10);
    plotDA(sim, figno,0);
    figno = figno + 1;
end

if sim.setup.plot_currents
    
    Jn = coeffs2sbp(sim.coeffs.Jn, sim);
    Jp = coeffs2sbp(sim.coeffs.Jp, sim);
    
    sim.setup.lp =0;
    
    makeplot(Jn, 'Currents', sim.phys.nmLs, sim.phys.mAicm2Js, sim, figno, 1,'red',10);
    
    makeplot(Jp, '', sim.phys.nmLs, sim.phys.mAicm2Js, sim, figno, 0,'blue',10);
    makeplot((Jn+Jp), '', sim.phys.nmLs, sim.phys.mAicm2Js, sim, figno, 0,'black',10);
    
    try
        for i = 1:(sim.setup.nx -1)
            Jn_hatr(i) = polyval(Jn(:,i),1) + sim.global.taunr(i) * (polyval(n(:,i), 1) - sim.coeffs.n_hat(i));
            Jn_hatl(i) = polyval(Jn(:,i+1),-1) + sim.global.taunl(i) * (polyval(n(:,i+1), -1) - sim.global.DeltaChi(i,i)*sim.coeffs.n_hat(i));
            
            Jp_hatr(i) = polyval(Jp(:,i),1) + sim.global.taupr(i) * (polyval(p(:,i), 1) - sim.coeffs.p_hat(i));
            Jp_hatl(i) = polyval(Jp(:,i+1),-1) + sim.global.taupl(i) * (polyval(p(:,i+1), -1) - sim.global.DeltaChiEg(i,i)*sim.coeffs.p_hat(i));
        end
        
        makeplot(Jn_hatr, '', sim.phys.nmLs, sim.phys.mAicm2Js, sim, figno, 0,'r.',0);
        makeplot(Jp_hatr, '', sim.phys.nmLs, sim.phys.mAicm2Js, sim, figno, 0,'b.',0);        
        makeplot(Jn_hatl, '', sim.phys.nmLs, sim.phys.mAicm2Js, sim, figno, 0,'r.',0);
        makeplot(Jp_hatl, '', sim.phys.nmLs, sim.phys.mAicm2Js, sim, figno, 0,'b.',0);
        
    catch
    end
    
    
    
    figno = figno + 1;
end

if sim.setup.plot_fields
    E = coeffs2sbp(sim.coeffs.E, sim);
    phi = coeffs2sbp(sim.coeffs.phi, sim);
    phin = phi + sim.sbp.dEg;
    phip = phin + sim.sbp.dChi;
    
    sim.setup.lp = 0;
    makeplot(E, 'E-field',  sim.phys.nmLs, 1, sim, figno, 1, 'green',10);
    makeplot(phi, 'Potential', sim.phys.nmLs, sim.phys.VVs, sim, figno, 0, 'cyan',10);
    makeplot(phin, 'Potential', sim.phys.nmLs, sim.phys.VVs, sim, figno, 0, 'b-',10);
    makeplot(phip, 'Potential', sim.phys.nmLs, sim.phys.VVs, sim, figno, 0, 'r-',10);
    makeplot(sim.coeffs.phi_hat', '', sim.phys.nmLs, sim.phys.VVs, sim, figno,0,'b.',0);
    
    
    figno = figno + 1;
end

if sim.setup.plot_levels
    plotlevels(sim, figno, 1)
    
    figno = figno + 1;
end

if sim.setup.plot_tau
    
    subplot(sim.setup.plotlayout(1), sim.setup.plotlayout(2), figno);
    
    sim.setup.lp = 0;
    try
        plot(sim.setup.x1, log10(sim.global.taunl), sim.setup.x1, log10(sim.global.taunr), sim.setup.x1, log10(sim.global.taupl), sim.setup.x1, log10(sim.global.taupr));
        legend('nl', 'nr','pl','pr');
        xlim([0,1]);
        ylim([-6,1.1]);
        %makeplot(1, 'tau_n',  sim.phys.nmLs, 1, sim, figno, 1, 'red',10);
        %makeplot(2, 'tau_p', sim.phys.nmLs, 1, sim, figno, 0, 'blue',10);
        %makeplot(3, 'tau_phi', sim.phys.nmLs, 1, sim, figno, 0, 'black',10);
    catch
    end
    
    figno = figno + 1;
end

if sim.setup.plot_U
    G = coeffs2sbp(sim.coeffs.G, sim);
    coeff = sbptimes(sim.sbp.ini2, sim.sbp.alpha, sim);
    n = coeffs2sbp(sim.coeffs.n, sim);
    p = coeffs2sbp(sim.coeffs.p, sim);
    np = sbptimes(n, p, sim);
    
%     n_sp = sbp2space(n, sim, 10);
%     n=n_sp(:,2);
%     p_sp = sbp2space(p, sim, 10);
%     p=p_sp(:,2);
%     ni_sp = sbp2space(sim.sbp.ni, sim, 10);
%     ni=ni_sp(:,2);
%     n1_sp = sbp2space(sim.sbp.n1, sim, 10);
%     n1=n1_sp(:,2)+0.01;
%     p1_sp = sbp2space(sim.sbp.p1, sim, 10);
%     p1=p1_sp(:,2);
%     taun_sp = sbp2space(sim.sbp.tausrhn, sim, 10);
%     taun=taun_sp(:,2);
%     taup_sp = sbp2space(sim.sbp.tausrhp, sim, 10);
%     taup=taup_sp(:,2);
%     
%     srhn = (n.*p-ni.*ni);
%     srhd = taup.*(p+p1)+taun.*(n+n1);
%     srh = srhn./srhd;
%     plot(srh);
    
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
    
    sim.setup.lp = 1;
    makeplot(G, 'G and R',  sim.phys.nmLs, sim.phys.icm3isGs, sim, figno, 1, 'green',10);
    makeplot(R, '', sim.phys.nmLs, sim.phys.icm3isGs, sim, figno, 0, 'red',10);
    
    makeplot(RD, '', sim.phys.nmLs, sim.phys.icm3isGs, sim, figno, 0, 'blue',10);
    makeplot(R_SRH, '', sim.phys.nmLs, sim.phys.icm3isGs, sim, figno, 0, 'yellow',10);
    
    makeplot(R_Aug, '', sim.phys.nmLs, sim.phys.icm3isGs, sim, figno, 0, 'black',10);
    
    makeplot(G - R, '', sim.phys.nmLs, sim.phys.icm3isGs, sim, figno, 0, 'red',10);
    figno = figno + 1;
end


drawnow;

%      makeplot(dphi, nmLs, cmLs/VVth, sim, 2, 1, 'black',10);
%      xlabel('z (nm)')
%      ylabel('phi (V)')



%Plot currents
%     Jndrift = polymatch(PN.mun,sim.pdeg).*n.*dphi;
%     Jndiff = Jn - Jndrift;
%
%     Jpdrift = polymatch(PN.mup,sim.pdeg).*p.*dphi;
%     Jpdiff = Jp - Jpdrift;
%
%
%     makeplot(Jndrift , nmLs, 1, sim, 4, 1, 'red',10);
%     makeplot(Jndiff , nmLs, 1, sim, 4, 0, 'r--',10);
%     makeplot(Jn , nmLs, 1, sim, 4, 0, 'r--',10);
%
%
%     makeplot(Jpdrift , nmLs, 1, sim, 4, 0, 'blue',10);
%     makeplot(Jpdiff , nmLs, 1, sim, 4, 0, 'b--',10);
%     makeplot(Jp, nmLs, 1, sim, 4, 0, 'b--',10);
%
%     makeplot(Jn+Jp, nmLs, 1, sim, 4, 1, 'black',10);
%
%     xlabel('z (nm)')
%     ylabel('phi (V)')

end