function [] = plotDA(sim, figno, clearplot)

%Depletion Approximation (valid if PN junction)
VVbi = sim.phys.VVs*sim.setup.phi_bi; % Vbi in eVs
ND = abs(polyval(sim.sbp.ND(:,1),-1));
NA = abs(polyval(sim.sbp.ND(:,end),1));

n0 = sim.phys.icm3Ns * sim.bc.n0;
n1 = sim.phys.icm3Ns * sim.bc.n1;
p0 = sim.phys.icm3Ns * sim.bc.p0;
p1 = sim.phys.icm3Ns * sim.bc.p1;

% Note that this assumes epsdc is constant
eps = sim.material.epsdc(1);
nmwscr = nm_from_cm(sqrt((2*eps*sim.phys.CiVicmeps0/sim.phys.Cq)*abs(VVbi)*(1/ND+1/NA)/sim.phys.icm3Ns));

nmwp = nmwscr * ND/(NA+ND);
nmwn = nmwscr * NA/(NA+ND);

nmLp = sim.material.nmSec(1);
nmLn = sim.material.nmSec(end);

nmz = -nmLp:1:nmLn;
Vphi_dep = 0 * nmz;
Vphi_dep(nmz<=-nmwp) = sim.phys.VVs*n1;
Vphi_dep(nmz>-nmwp) = sim.phys.VVs*n1 + 1e-14 * sim.phys.Cq*sim.phys.icm3Ns*NA/(2*eps*sim.phys.CiVicmeps0)*(nmz(nmz>-nmwp)+nmwp).^2; % 1e-14 converts z and wp to cm to cancel NA in cm^-2
Vphi_dep(nmz>0) = sim.phys.VVs*n1 + VVbi - 1e-14 * sim.phys.Cq*sim.phys.icm3Ns*ND/(2*eps*sim.phys.CiVicmeps0)*(nmz(nmz>0)-nmwn).^2; %
Vphi_dep(nmz>=nmwn) = sim.phys.VVs*n1 + VVbi;


n_dep = 0 * nmz;
n_dep(nmz<=-nmwp) = n0;
n_dep(nmz>-nmwp) = n1; % 1e-14 converts z and wp to cm to cancel NA in cm^-2
n_dep(nmz>0) = n1; %
n_dep(nmz>=nmwn) = n1;

p_dep = 0 * nmz;
p_dep(nmz<=-nmwp) = p0;
p_dep(nmz>-nmwp) = p0; % 1e-14 converts z and wp to cm to cancel NA in cm^-2
p_dep(nmz>0) = p0; %
p_dep(nmz>=nmwn) = p1;

subplot(sim.setup.plotlayout(1), sim.setup.plotlayout(2), figno);

if ~clearplot
    hold on
end
if sim.setup.lp
    plot(nmz + nmLp, log10(n_dep) ,'red--') % Space field
    plot(nmz + nmLp, log10(p_dep) ,'blue--') % Space field
else
        plot(nmz + nmLp, n_dep ,'red--') % Space field
    plot(nmz + nmLp, p_dep ,'blue--') % Space field
end
    hold off
            xlabel('z (nm)')
    ylabel('n,p,ni (cm^{-3})')

end
