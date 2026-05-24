function [] = plotlevels(sim, plotno, clearplot) 

    chidata =  sbp2space(sim.sbp.Chi, sim)';
    xdata = sim.phys.nmLs * chidata(1,:);
        
    Evacdata = -coeffs2space(sim.coeffs.phi, sim)';
    ndata = coeffs2space(sim.coeffs.n, sim)';
    pdata = coeffs2space(sim.coeffs.p, sim)';
    nidata = sbp2space(sim.sbp.ni, sim)';
    Egdata = sbp2space(sim.sbp.Eg, sim)';

    Ncdata = sbp2space(sim.sbp.Nc,sim)';
    Nvdata = sbp2space(sim.sbp.Nv,sim)';
    
    
    % Distance from E_fn to phi at LHS
    Vc0 = - log(sim.bc.n0/Ncdata(2,1));
    EF0 = Evacdata(2,1) - chidata(2,1) - Vc0;
    
    Evac1 = [Evacdata(1,:); Evacdata(2,:)];
    Ec = [Evacdata(1,:); Evacdata(2,:) - chidata(2,:)];
    Ei = [Evacdata(1,:); Ec(2,:) - 0.5*Egdata(2,:) - 0.5*log(Ncdata(2,:)./Nvdata(2,:))];
    Fn = [Evacdata(1,:); Ei(2,:) + log(ndata(2,:)./nidata(2,:))];
    Fp = [Evacdata(1,:); Ei(2,:) - log(pdata(2,:)./nidata(2,:))];  
    Ev = [Evacdata(1,:); Ec(2,:) - Egdata(2,:)];
    
            
    
    subplot(sim.setup.plotlayout(1), sim.setup.plotlayout(2), plotno);
    
    if clearplot
        hold off
    else
        hold on
    end
    plot(xdata(1,:), sim.phys.VVs * Evac1(2,:),'black') % Space field
    
    hold on
    plot(xdata(1,:), sim.phys.VVs * real(Ec(2,:)) , 'red') % Conduction Band
    plot(xdata(1,:), sim.phys.VVs * real(Fn(2,:)), 'red--') % Electron Fermi level
    plot(xdata(1,:), sim.phys.VVs * real(Ei(2,:)), 'black--') % Intrinsic energy level
    plot(xdata(1,:), sim.phys.VVs * real(Fp(2,:)), 'blue--') % Hole Fermi level
    plot(xdata(1,:), sim.phys.VVs * real(Ev(2,:)), 'blue') % Valence Band
 
    xlabel('z (nm)')
    ylabel('Potentials (V)')
    
    legend('Phi', 'E_c', 'F_n', 'F_p', 'E_v');
    hold off
    
end