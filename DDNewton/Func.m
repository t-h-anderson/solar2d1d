

function [ sim ] = Func(sim, varargin)

% Calculate Tensor products
if sim.setup.tensor
    sim.global.Psi_phi = tiprod(sim.global.Psi, 3, sim.coeffs.phi).'; % Gives A_ji
    sim.global.Psi_n = tiprod(sim.global.Psi, 1, sim.coeffs.n); % Gives A_jk
    sim.global.Psi_p = tiprod(sim.global.Psi, 1, sim.coeffs.p); % Gives A_jk
else
    sim.global.Psi_phi = tiprods(sim.global.Psis3, sim.coeffs.phi, sim).'; % Gives A_ji
    sim.global.Psi_n = tiprods(sim.global.Psis1, sim.coeffs.n, sim).'; % Gives A_ji
    sim.global.Psi_p = tiprods(sim.global.Psis1, sim.coeffs.p, sim).'; % Gives A_ji
    
    sim.global.Psialpha_n = tiprods(sim.global.Psialphas1, sim.coeffs.n, sim).'; % Gives A_ji
    sim.global.Psialpha_p = tiprods(sim.global.Psialphas3, sim.coeffs.p, sim).'; % Gives A_ji
end

% Calculate recombination terms using quadrature of mass lumping
if sim.setup.quadtoggle == 1

    % Quadrature integration of Auger and SRH terms
    sim.global.R_Aug = sim.setup.lt4coeffs .* QuadIntVect(@(x) x(:,4) .* x(:,1) .* (x(:,1) .* x(:,2) - x(:,3).*x(:,3)) + x(:,5) .* x(:,2) .* (x(:,1) .* x(:,2) - x(:,3).*x(:,3)),...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.Cn, sim.coeffs.Cp], sim);
sim.global.dR_Augdn = sim.setup.lt4coeffs .* QuadIntMatrix(@(x) x(:,4) .* (2 * x(:,1) .* x(:,2) - x(:,3).*x(:,3)) + x(:,5) .* x(:,2) .^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.Cn, sim.coeffs.Cp], sim);
sim.global.dR_Augdp = sim.setup.lt4coeffs .* QuadIntMatrix(@(x) x(:,5) .* (2 * x(:,1) .* x(:,2) - x(:,3).*x(:,3)) + x(:,4) .* x(:,1) .^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.Cn, sim.coeffs.Cp], sim);
    
sim.global.R_SRH = sim.setup.damping * sim.setup.lt4coeffs .* QuadIntVect(@(x)(x(:,1) .* x(:,2) - x(:,3).*x(:,3))./(x(:,5).*(x(:,1) + x(:,6)) + x(:,4).*(x(:,2) + x(:,7))),...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.tausrhn, sim.coeffs.tausrhp, sim.coeffs.n1, sim.coeffs.p1], sim);
sim.global.dR_SRHdn = sim.setup.damping * sim.setup.lt4coeffs .* QuadIntMatrix(@(x)(x(:,4) .* x(:,2) .* (x(:,2)+x(:,7)) + x(:,5) .* (x(:,2) .* x(:,6) +x(:,3).*x(:,3))) ./ (x(:,5) .* (x(:,1) + x(:,6)) + x(:,4) .* (x(:,2) + x(:,7))).^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.tausrhn, sim.coeffs.tausrhp, sim.coeffs.n1, sim.coeffs.p1], sim);
sim.global.dR_SRHdp = sim.setup.damping * sim.setup.lt4coeffs .* QuadIntMatrix(@(x)(x(:,5) .* x(:,1) .* (x(:,1)+x(:,6)) + x(:,4) .* (x(:,1) .* x(:,7) +x(:,3).*x(:,3))) ./ (x(:,5) .* (x(:,1) + x(:,6)) + x(:,4) .* (x(:,2) + x(:,7))).^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.tausrhn, sim.coeffs.tausrhp, sim.coeffs.n1, sim.coeffs.p1], sim);

else
    % Mass lumping integration of Auger and SRH terms
    sim.global.R_Aug = sim.setup.lt4coeffs .* LobInt(@(x) x(:,4) .* x(:,1) .* (x(:,1) .* x(:,2) - x(:,3).*x(:,3)) + x(:,5) .* x(:,2) .* (x(:,1) .* x(:,2) - x(:,3).*x(:,3)),...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.Cn, sim.coeffs.Cp], sim);
sim.global.dR_Augdn = diag(sim.setup.lt4coeffs .* LobInt(@(x) x(:,4) .* (2 * x(:,1) .* x(:,2) - x(:,3).*x(:,3)) + x(:,5) .* x(:,2) .^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.Cn, sim.coeffs.Cp], sim));
sim.global.dR_Augdp = diag(sim.setup.lt4coeffs .* LobInt(@(x) x(:,5) .* (2 * x(:,1) .* x(:,2) - x(:,3).*x(:,3)) + x(:,4) .* x(:,1) .^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.Cn, sim.coeffs.Cp], sim));

    
sim.global.R_SRH = sim.setup.damping * sim.setup.lt4coeffs .* LobInt(@(x)(x(:,1) .* x(:,2) + x(:,3).*x(:,3))./(x(:,5).*(x(:,1) + x(:,6)) + x(:,4).*(x(:,2) + x(:,7))),...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.tausrhn, sim.coeffs.tausrhp, sim.coeffs.n1, sim.coeffs.p1], sim);
sim.global.dR_SRHdn = diag(sim.setup.damping * sim.setup.lt4coeffs .* LobInt(@(x)(x(:,4) .* x(:,2) .* (x(:,2)+x(:,7)) + x(:,5) .* (x(:,2) .* x(:,6) +x(:,3).*x(:,3))) ./ (x(:,5) .* (x(:,1) + x(:,6)) + x(:,4) .* (x(:,2) + x(:,7))).^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.tausrhn, sim.coeffs.tausrhp, sim.coeffs.n1, sim.coeffs.p1], sim));
sim.global.dR_SRHdp = diag(sim.setup.damping * sim.setup.lt4coeffs .* LobInt(@(x)(x(:,5) .* x(:,1) .* (x(:,1)+x(:,6)) + x(:,4) .* (x(:,1) .* x(:,7) +x(:,3).*x(:,3))) ./ (x(:,5) .* (x(:,1) + x(:,6)) + x(:,4) .* (x(:,2) + x(:,7))).^2,...
    [sim.coeffs.n, sim.coeffs.p, sim.coeffs.ni, sim.coeffs.tausrhn, sim.coeffs.tausrhp, sim.coeffs.n1, sim.coeffs.p1], sim));
end

% Recombination and total generation terms
R = (sim.global.Psialpha_n*sim.coeffs.p - sim.global.alpha) + sim.global.R_SRH + sim.global.R_Aug;
U = sim.setup.damping * (sim.global.G - R);

% Build RHS
Fn1 = sim.global.Pn * sim.coeffs.Jn + (sim.global.Psi_phi + sim.global.Vn + sim.global.P01) * sim.coeffs.n + sim.global.S_psi_n.' * sim.coeffs.n_hat;
Fn2 = (sim.global.D1-sim.global.P01) * sim.coeffs.Jn + sim.global.D_taun * sim.coeffs.n + U + sim.global.S_taunpsi_n.' * sim.coeffs.n_hat;
Fn_hat = sim.global.S_psi * sim.coeffs.Jn + sim.global.S_taunpsi * sim.coeffs.n - sim.global.S_taun_n * sim.coeffs.n_hat;

Fp1 = sim.global.Pp * sim.coeffs.Jp + (- sim.global.Psi_phi - sim.global.Vp + sim.global.P01) * sim.coeffs.p + sim.global.S_psi_p.' *sim.coeffs.p_hat;
Fp2 = (sim.global.D1-sim.global.P01) * sim.coeffs.Jp + sim.global.D_taup * sim.coeffs.p - U + sim.global.S_tauppsi_p.' * sim.coeffs.p_hat;
Fp_hat = sim.global.S_psi * sim.coeffs.Jp + sim.global.S_tauppsi*sim.coeffs.p - sim.global.S_taup_p'*sim.coeffs.p_hat;

Fphi1 = sim.global.Pphi * sim.coeffs.E + sim.global.P01 * sim.coeffs.phi + sim.global.S_psi.' * sim.coeffs.phi_hat;
Fphi2 = (sim.global.D1-sim.global.P01) * sim.coeffs.E + sim.global.P00 * (sim.coeffs.n - sim.coeffs.p - sim.coeffs.ND) + sim.global.D_tauphi * sim.coeffs.phi  + sim.global.S_tauphipsi.' * sim.coeffs.phi_hat;
Fphi_hat = sim.global.S_psi * sim.coeffs.E + sim.global.S_tauphipsi * sim.coeffs.phi - sim.global.S_tauphi.' * sim.coeffs.phi_hat;

F = [Fn1; Fn2; Fn_hat; Fp1; Fp2; Fp_hat; Fphi1; Fphi2; Fphi_hat];

% Add boundary conditions
sim.global.F = F + sim.global.B;
end



