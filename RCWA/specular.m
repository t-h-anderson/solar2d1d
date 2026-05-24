
function [absorb]=specular...
    (ll, tt, sim)

%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computations of transmission \
% coefficients, TP[i] contains the tranmission coefficients at the ith
% interface so TP[Ns+1] contains the coefficients of whole structure as
% their are Ns+1 interfaces

inmk0 = sim.phys.inmk0(ll);
Nt = sim.setup.Nt;
epsf = sim.zfl.eps{ll};
radtheta = sim.phys.radtheta(tt);
n0 = sqrt(epsf(1, 2* sim.setup.Nt +1));


inmkxn = n0 * inmk0.*sin(radtheta)+(-Nt:Nt)*2*pi/sim.input.nmLx;
inmkzn = sqrt(n0^2.*inmk0^2-inmkxn.^2);

RP2 = abs(sim.results.R).^2;    % Reflection amplitudes
TP2 = abs(sim.results.T).^2;    % Transmission amplitudes


% The following block computes a diagonal matrix with real parts of \
% kz^(n) as its elements**)
kzr1 = real(inmkzn)/(n0*inmk0*cos(radtheta));
kzr2 = real(inmkzn)/(n0*inmk0*cos(radtheta));
RKZ1 = diag(kzr1);
RKZ2 = diag(kzr2);

%  (******************************* p incident, p reflected and \
% transmitted Pp ***************************************************

RPp = (RP2'*RKZ1).';
TPp = (TP2'*RKZ2).';
rp0 = RPp(Nt+1);
rp1 = sum(RPp)-RPp(Nt+1);
tp0 = TPp(Nt+1);
tp1 = sum(TPp)-TPp(Nt+1);
% (***************************************************************

absorb = 1-rp0-rp1-tp0-tp1;
end
