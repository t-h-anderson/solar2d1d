
function [sim] = setboundaries(sim)

% Set boundary condition vector B
np = sim.setup.np;
pdeg1 = sim.setup.pdeg1;
nx = sim.setup.nx;
sim.global.B = sparse(6*np+3*(sim.setup.nx-1),1);

for i=1:pdeg1
    sim.global.B(i) = sim.global.B(i) + sim.bc.n0*sim.local.pl(i);
    j = np - pdeg1;
    sim.global.B(j+i) = sim.global.B(j+i) - sim.bc.n1*sim.local.pr(i);

    j = np;
    sim.global.B(i+j) = sim.global.B(i+j) + sim.global.taun0*sim.bc.n0*sim.local.pl(i);
    j = 2*np - pdeg1;
    sim.global.B(i+j) = sim.global.B(i+j) - sim.global.taun1*sim.bc.n1*sim.local.pr(i);

    j = 2*np+(nx-1);
    sim.global.B(i+j) = sim.global.B(i+j) + sim.bc.p0*sim.local.pl(i);
    j = 3*np+(nx-1) - pdeg1;
    sim.global.B(i+j) = sim.global.B(i+j) -  sim.bc.p1*sim.local.pr(i);

    j = 3*np+(nx-1);
    sim.global.B(i+j) = sim.global.B(i+j) + sim.global.taup0*sim.bc.p0*sim.local.pl(i);
    j = 4*np+(nx-1) - pdeg1;
    sim.global.B(i+j) = sim.global.B(i+j) - sim.global.taup1*sim.bc.p1*sim.local.pr(i);
    
    j = 4*np+2*(nx-1);
    sim.global.B(i+j) = sim.global.B(i+j) + sim.bc.phi0*sim.local.pl(i);
    j = 5*np+2*(nx-1) - pdeg1;
    sim.global.B(i+j) = sim.global.B(i+j) - sim.bc.phi1*sim.local.pr(i);

    j = 5*np+2*(nx-1);
    sim.global.B(i+j) = sim.global.B(i+j) + sim.global.tauphi0*sim.bc.phi0*sim.local.pl(i);
    j = 6*np+2*(nx-1) - pdeg1;
    sim.global.B(i+j) = sim.global.B(i+j) - sim.global.tauphi1*sim.bc.phi1*sim.local.pr(i);
    
end



pos = 0;
sim.global.Bn1 = sim.global.B(pos + 1:pos + np);

pos = pos + np;
sim.global.Bn2 = sim.global.B(pos + 1: pos + np);

pos = pos + np;
sim.global.Bn_hat = sim.global.B(pos + 1: pos + nx - 1);

pos = pos + nx - 1;
sim.global.Bp1 = sim.global.B(pos + 1:pos + np);

pos = pos + np;
sim.global.Bp2 = sim.global.B(pos + 1 : pos + np);

pos = pos + np;
sim.global.Bp_hat = sim.global.B(pos + 1: pos + nx - 1);

pos = pos + nx - 1;
sim.global.Bphi1 = sim.global.B(pos + 1:pos + np);

pos = pos + np;
sim.global.Bphi2 = sim.global.B(pos + 1 : pos + np);

pos = pos + np;
sim.global.Bphi_hat = sim.global.B(pos + 1: pos + nx - 1);

end