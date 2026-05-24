function [ sim ] = BCond(sim)
% Set the initial conditions for the simulation

n_te = zeros(sim.setup.pdeg1, 2);
p_te = zeros(sim.setup.pdeg1, 2);

for i =1:2
    
    if i == 1
        ND = sim.sbp.ND(:,1);
        ratio = polyval(sim.sbp.ni(:,1),0)/(polyval(ND,0));
        ni2 = sim.sbp.ni2(:,1);
        ND2 = sim.sbp.ND2(:,1);
        iND = sim.sbp.iND(:,1);       
    else
        ND =  sim.sbp.ND(:,end);
        ratio = polyval(sim.sbp.ni(:,end),0)/(polyval(ND,0));
        ni2 = sim.sbp.ni2(:,end);
        ND2 = sim.sbp.ND2(:,end);
        iND = sim.sbp.iND(:,end);
    end
    
    
    if  (ratio) < 1e-4 && (ratio) >=0
        n_te(:,i) =  ND;
        p_te(:,i) = polymatch(polytimes(ni2, iND), sim.setup.pdeg1);
    elseif ratio > -1e-4 && ratio <=0
        p_te(:,i) = -ND;
        n_te(:,i) = - polymatch(polytimes(ni2, iND), sim.setup.pdeg1);
    else
        n_te(:,i) = (ND +  polyfunc(@(x) sqrt(x), (ND2 + 4*ni2)))/2;
        p_te(:,i) = (-ND +  polyfunc(@(x) sqrt(x), (ND2 + 4*ni2)))/2;
    end
    
end


% Calculate Dirichelt boundary conditions for n
sim.bc.n0 = polyval(n_te(:, 1), -1);
sim.bc.n1 = polyval(n_te(:, end), 1);

% Calculate Dirichelt boundary conditions for p
sim.bc.p0 = polyval(p_te(:, 1), -1);
sim.bc.p1 = polyval(p_te(:, end), 1);

if imag(sim.bc.p0)~=0 || imag(sim.bc.p1)~=0 || imag(sim.bc.n0)~=0 || imag(sim.bc.n1)~=0
    pause(0)
end

% Calculate initial value for phi
Nc0 = polyval(sim.sbp.Nc(:,1), -1);
Chi0 = polyval(sim.sbp.Chi(:,1), -1);
Nc1 = polyval(sim.sbp.Nc(:,end), 1);
Chi1 = polyval(sim.sbp.Chi(:,end), 1);

sim.bc.phi0 = 0;
sim.bc.phi1 = ((Chi0 - Chi1) + (log(Nc0/sim.bc.n0) - log(Nc1/sim.bc.n1))) + sim.setup.Vext; % equivalent to version in paper

if imag(sim.bc.phi0)~=0 || imag(sim.bc.phi1)~=0
    pause(0);
end


end

