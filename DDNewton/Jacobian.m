
function [ J, sim ] = Jacobian( sim, delta )
% Calculates the Jacobian of our function

nx = sim.setup.nx;
pdeg = sim.setup.pdeg1;
np = nx*pdeg;



if delta == 0

    zeroS = sparse(nx-1, nx-1);
    zeroM1 = sparse(nx-1, np);
    zeroM2 = sparse(np, nx-1);
    zeroL = sparse(np, np);
    
    dRn = sim.global.Psialpha_n + sim.global.dR_SRHdn + sim.global.dR_Augdn;
    dRp = sim.global.Psialpha_p + sim.global.dR_SRHdp + sim.global.dR_Augdp;
    
    J = [...
        [sim.global.Psi_phi + sim.global.Vn + sim.global.P01; sim.global.D_taun - dRn ; sim.global.S_taunpsi; zeroL; dRn; zeroM1; zeroL; sim.global.P00; zeroM1],...
        [sim.global.Pn; sim.global.D1 - sim.global.P01; sim.global.S_psi; zeroL; zeroL; zeroM1; zeroL; zeroL; zeroM1],...
        [ sim.global.S_psi_n.'; sim.global.S_taunpsi_n.'; -sim.global.S_taun_n; zeroM2; zeroM2; zeroS; zeroM2; zeroM2; zeroS] ,...
        [zeroL; -dRp; zeroM1; - sim.global.Psi_phi - sim.global.Vp + sim.global.P01; sim.global.D_taup + dRp ; sim.global.S_tauppsi; zeroL; - sim.global.P00; zeroM1],...
        [zeroL; zeroL; zeroM1; sim.global.Pp; sim.global.D1 - sim.global.P01; sim.global.S_psi; zeroL; zeroL; zeroM1],...
        [zeroM2; zeroM2; zeroS; sim.global.S_psi_p.'; sim.global.S_tauppsi_p.'; -sim.global.S_taup_p;  zeroM2; zeroM2; zeroS],...
        [sim.global.Psi_n; zeroL; zeroM1; - sim.global.Psi_p ; zeroL; zeroM1; sim.global.P01; sim.global.D_tauphi; sim.global.S_tauphipsi],...
        [zeroL; zeroL; zeroM1; zeroL; zeroL; zeroM1; sim.global.Pphi; sim.global.D1 - sim.global.P01; sim.global.S_psi],...
        [zeroM2; zeroM2; zeroS;  zeroM2; zeroM2; zeroS;  sim.global.S_psi.';  sim.global.S_tauphipsi.'; -sim.global.S_tauphi]...
        ];
    
elseif delta > 0
    
    J = zeros(6*np + 3*(nx-1));
    for i = 1:np
        
        a = 0; b=0;
        shift = a*np + b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.n(i) = sim2.coeffs.n(i)+delta;
        sim2 = Func(sim2);
        F2n = sim2.global.F;
        dFn = (F2n - sim.global.F)/delta;
        
        J(:, i + shift) = dFn;
        
        a = 1; b=0;
        shift = a * np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.Jn(i) = sim2.coeffs.Jn(i)+delta;
        sim2 = Func(sim2);
        F2Jn = sim2.global.F;
        dFJn = (F2Jn -sim.global.F)/delta;
        
        J(:, i + shift) = dFJn;
        
        a = 2; b=1;
        shift = a*np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.p(i) = sim2.coeffs.p(i)+delta;
        sim2 = Func(sim2);
        F2p = sim2.global.F;
        dFp = (F2p -sim.global.F)/delta;
        
        J(:, i + shift) = dFp;
        
        a = 3; b=1;
        shift = a*np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.Jp(i) = sim2.coeffs.Jp(i)+delta;
        sim2 = Func(sim2);
        F2Jp = sim2.global.F;
        dFJp = (F2Jp -sim.global.F)/delta;
        
        J(:, i + shift) = dFJp;
        
        a = 4; b=2;
        shift = a*np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.phi(i) = sim2.coeffs.phi(i)+delta;
        sim2 = Func(sim2);
        F2phi = sim2.global.F;
        dFphi = (F2phi -sim.global.F)/delta;
        
        J(:, i + shift) = dFphi;
        
        a = 5; b=2;
        shift = a*np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.E(i) = sim2.coeffs.E(i)+delta;
        sim2 = Func(sim2);
        F2E = sim2.global.F;
        dFE = (F2E -sim.global.F)/delta;
        
        J(:, i + shift) = dFE;
    end
    
    
    for i = 1:(nx-1)
        
        a = 2; b=0;
        shift = a*np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.n_hat(i) = sim2.coeffs.n_hat(i)+delta;
        sim2 = Func(sim2);
        F2n_hat = sim2.global.F;
        dFn_hat = (F2n_hat -sim.global.F)/delta;
        
        J(:, i + shift) = dFn_hat;
        
        a = 4; b=1;
        shift = a*np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.p_hat(i) = sim2.coeffs.p_hat(i)+delta;
        sim2 = Func(sim2);
        F2p_hat = sim2.global.F;
        dFp_hat = (F2p_hat -sim.global.F)/delta;
        
        J(:, i + shift) = dFp_hat;
        
        a = 6; b=2;
        shift = a*np+ b*(nx - 1);
        sim2 = sim;
        sim2.coeffs.phi_hat(i) = sim2.coeffs.phi_hat(i)+delta;
        sim2 = Func(sim2);
        F2phi_hat = sim2.global.F;
        dFphi_hat = (F2phi_hat -sim.global.F)/delta;
        
        J(:, i + shift) = dFphi_hat;
        
        
    end
    
else
    
    error('Negative delta unused');
    zeroS = sparse(nx-1, nx-1);
    zeroM1 = sparse(nx-1, np);
    zeroM2 = sparse(np, nx-1);
    zeroL = sparse(np, np);
    
    GRdamping = 0;
    
    J = [[sim.global.Psi_phi + sim.global.Vn + sim.global.P01; sim.global.D_taun; sim.global.S_taunpsi; zeroL; zeroL; zeroM1; zeroL; sim.global.P00; zeroM1],...
        [sim.global.Pn - GRdamping * sim.global.Psialpha_p; sim.global.D1 - sim.global.P01; sim.global.S_psi; zeroL; zeroL; zeroM1; zeroL; zeroL; zeroM1],...
        [ sim.global.S_psi.'; sim.global.S_taunpsi.'; -sim.global.S_taun; zeroM2; zeroM2; zeroS; zeroM2; zeroM2; zeroS] ,...
        [zeroL; zeroL; zeroM1; - sim.global.Psi_phi - sim.global.Vp + sim.global.P01; sim.global.D_taup; sim.global.S_tauppsi; zeroL; - sim.global.P00; zeroM1],...
        [zeroL; zeroL; zeroM1; sim.global.Pp + GRdamping * sim.global.Psialpha_n; sim.global.D1 - sim.global.P01; sim.global.S_psi; zeroL; zeroL; zeroM1],...
        [zeroM2; zeroM2; zeroS; sim.global.S_psi.'; sim.global.S_tauppsi.'; -sim.global.S_taup;  zeroM2; zeroM2; zeroS],...
        [sim.global.Psi_n; zeroL; zeroM1; - sim.global.Psi_p ; zeroL; zeroM1; sim.global.P01; sim.global.D_tauphi; sim.global.S_tauphipsi],... % Check
        [zeroL; zeroL; zeroM1; zeroL; zeroL; zeroM1; sim.global.Pphi; sim.global.D1 - sim.global.P01; sim.global.S_psi],...
        [zeroM2; zeroM2; zeroS;  zeroM2; zeroM2; zeroS;  sim.global.S_psi.';  sim.global.S_tauphipsi.'; -sim.global.S_tauphi ]];
    
end
end




