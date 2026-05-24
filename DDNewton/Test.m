
function sim = Test(sim)

delta = 1e-7;

nx = sim.setup.nx;
pdeg1 = sim.setup.pdeg1;
np = sim.setup.np;

for i = 1:(nx*pdeg1)

    a = 0; b=0;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.n(i) = sim2.coeffs.n(i)+delta;
    sim2 = Func(sim2);
    F2n = sim2.global.F;
    dFn = (F2n - sim.global.F)/delta;
    dFne = sim.global.J(:, i + shift);
    diffn = dFn - dFne;
    errorn(i) = max(abs(diffn));
    
    a = 1; b=0;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.Jn(i) = sim2.coeffs.Jn(i)+delta;
    sim2 = Func(sim2);
    F2Jn = sim2.global.F;
    dFJn = (F2Jn -sim.global.F)/delta;
    diffJn = dFJn -sim.global.J(:,i + shift);
    errorJn(i) = max(abs(diffJn));

    
     a = 2; b=1;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.p(i) = sim2.coeffs.p(i)+delta;
    sim2 = Func(sim2);
    F2p = sim2.global.F;
    dFp = (F2p -sim.global.F)/delta;
    dFpe = sim.global.J(:, i + shift);
    diffp = dFp - dFpe;
    errorp(i) = max(abs(diffp));
    
    a = 3; b=1;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.Jp(i) = sim2.coeffs.Jp(i)+delta;   
    sim2 = Func(sim2);
    F2Jp = sim2.global.F;
    dFJp = (F2Jp -sim.global.F)/delta;
    diffJp = dFJp -sim.global.J(:,i + shift);
    errorJp(i) = max(abs(diffJp));

     a = 4; b=2;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.phi(i) = sim2.coeffs.phi(i)+delta;
    sim2 = Func(sim2);
    F2phi = sim2.global.F;
    dFphi = (F2phi -sim.global.F)/delta;
    diffphi = dFphi -sim.global.J(:,i + shift);
    errorphi(i) = max(abs(diffphi));
    
    a = 5; b=2;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.E(i) = sim2.coeffs.E(i)+delta;
    sim2 = Func(sim2);
    F2E = sim2.global.F;
    dFE = (F2E -sim.global.F)/delta;
    diffE = dFE -sim.global.J(:,i + shift);
    errorE(i) = max(abs(diffE));
end


for i = 1:(nx-1)

    a = 2; b=0;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.n_hat(i) = sim2.coeffs.n_hat(i)+delta;
    sim2 = Func(sim2);
    F2n_hat = sim2.global.F;
    dFn_hat = (F2n_hat -sim.global.F)/delta;
    diffn_hat = dFn_hat -sim.global.J(:,i + shift);
    errorn_hat(i) = max(abs(diffn_hat));
    
    
     a = 4; b=1;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.p_hat(i) = sim2.coeffs.p_hat(i)+delta;
    sim2 = Func(sim2);
    F2p_hat = sim2.global.F;
    dFp_hat = (F2p_hat -sim.global.F)/delta;
    diffp_hat = dFp_hat -sim.global.J(:,i + shift);
    errorp_hat(i) = max(abs(diffp_hat));
    
    a = 6; b=2;
    shift = a*np + b*(nx - 1);    
    sim2 = sim;
    sim2.coeffs.phi_hat(i) = sim2.coeffs.phi_hat(i)+delta;
    sim2 = Func(sim2);
    F2phi_hat = sim2.global.F;
    dFphi_hat = (F2phi_hat -sim.global.F)/delta;
    diffphi_hat = dFphi_hat -sim.global.J(:,i + shift);
    errorphi_hat(i) = max(abs(diffphi_hat));
        
end


errorm = max(abs([errorn, errorp, errorphi]));
errorJ = max(abs([errorJn, errorJp, errorE]));
error_hat = max(abs([errorn_hat, errorp_hat, errorphi_hat]));

errormax = abs(max([errorm,errorJ,error_hat]));

if isfield(sim.results,'errormax')
    sim.results.errormax = [sim.results.errormax, errormax];
else
    sim.results.errormax = errormax;
end
if errormax > 1e-3
   display('Error > 1e-3');
   %error('Jacobian test higher than tolerance');
end

end


