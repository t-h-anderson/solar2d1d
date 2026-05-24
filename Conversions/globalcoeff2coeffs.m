function sim = globalcoeff2coeffs(xnew, sim)
 
np = sim.setup.np;
nx = sim.setup.nx;

pos = 0;
sim.coeffs.n = xnew(pos + 1:pos + np);

pos = pos + np;
sim.coeffs.Jn = xnew(pos + 1: pos + np);

pos = pos + np;
sim.coeffs.n_hat = xnew(pos + 1: pos + nx - 1);

pos = pos + nx - 1;
sim.coeffs.p = xnew(pos + 1:pos + np);

pos = pos + np;
sim.coeffs.Jp = xnew(pos + 1 : pos + np);

pos = pos + np;
sim.coeffs.p_hat = xnew(pos + 1: pos + nx - 1);

pos = pos + nx - 1;
sim.coeffs.phi = xnew(pos + 1:pos + np);

pos = pos + np;
sim.coeffs.E = xnew(pos + 1 : pos + np);

pos = pos + np;
sim.coeffs.phi_hat = xnew(pos + 1: pos + nx - 1);

end