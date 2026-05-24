
function [ sim ] = initializeSim( sim )
% sim = initialiseSim(sim)
% Initialises all of the permatnent matrices in the simulation.
% This is independent of the junction design.

if sim.setup.SolverSteps
    disp('Initialising');
end

sim.setup.damping = sim.setup.dampmin;

% Structure the output
nplot = sim.setup.plot_densities + sim.setup.plot_fields + sim.setup.plot_levels +sim.setup.plot_currents + sim.setup.plot_U + sim.setup.plot_tau;

ploty = max(ceil(sqrt(nplot)),1);
plotx = max(ceil(nplot/ploty),1);

sim.setup.plotlayout = [plotx, ploty];

% Create matrix for scalling integration over standard basis.
sim.setup.lt = (sim.setup.mesh)/2;
sim.setup.lt4coeffs =  kron(sim.setup.lt,ones(sim.setup.pdeg1,1));

% Create the mesh endpoints for [0,Lx]
sim.setup.x_r = tril(ones(sim.setup.nx))*sim.setup.mesh;
sim.setup.x_l = [0 ; sim.setup.x_r(1:end-1)];

% Create the coeficient variables
np = sim.setup.np;

sim.coeffs.n = zeros(np,1);
sim.coeffs.Jn = zeros(np,1);
sim.coeffs.p = zeros(np,1);
sim.coeffs.Jp = zeros(np,1);
sim.coeffs.phi = zeros(np,1);
sim.coeffs.E = zeros(np,1);

sim.coeffs.n_hat = zeros(sim.setup.nx-1,1);
sim.coeffs.p_hat = zeros(sim.setup.nx-1,1);
sim.coeffs.phi_hat = zeros(sim.setup.nx-1,1);

% Create matrixes to speed up integration
% These will need changeing if differnet poly sizes used for n, p, phi
pdeg1 = sim.setup.pdeg1;
qdeg1 = sim.setup.pdeg1;
sim.int.nnmat = zeros(pdeg1,qdeg1);
vect = zeros(2*max(pdeg1,qdeg1) - 1,1);
vect(1:2:2*max(pdeg1,qdeg1)-1) = 1./(2*max(pdeg1,qdeg1)-1:-2:0);
for i=1:qdeg1
    sim.int.nnmat(:,i) = vect(end - qdeg1 + i - pdeg1 + 1:end - qdeg1 + i);
end

pdeg1 = sim.setup.pdeg1;
qdeg1 = 2*sim.setup.pdeg+1;
sim.int.n2nmat = zeros(pdeg1,qdeg1);
vect = zeros(2*max(pdeg1,qdeg1) - 1,1);
vect(1:2:2*max(pdeg1,qdeg1)-1) = 1./(2*max(pdeg1,qdeg1)-1:-2:0);
for i=1:qdeg1
    sim.int.n2nmat(:,i) = vect(end - qdeg1 + i - pdeg1 + 1:end - qdeg1 + i);
end

pdeg1 = 2*sim.setup.pdeg+1;
qdeg1 = 2*sim.setup.pdeg+1;
sim.int.n2n2mat = zeros(pdeg1,qdeg1);
vect = zeros(2*max(pdeg1,qdeg1) - 1,1);
vect(1:2:2*max(pdeg1,qdeg1)-1) = 1./(2*max(pdeg1,qdeg1)-1:-2:0);
for i=1:qdeg1
    sim.int.n2n2mat(:,i) = vect(end - qdeg1 + i - pdeg1 + 1:end - qdeg1 + i);
end

% Integration weights for Lobatto integration
switch sim.setup.pdeg1
    case 3
        wi = [1/3, 4/3, 1/3].';
    
    case 4
        wi = [1/6, 5/6, 5/6, 1/6].';
        
    case 5
        wi = [1/10, 49/90, 32/45, 49/90, 1/10].';
        
    case 6
        wi = [1/15, (14-sqrt(7))/30, (14+sqrt(7))/30, (14+sqrt(7))/30, (14-sqrt(7))/30, 1/15].';
    otherwise
         [~,wi]=lobatto_compute(sim.setup.pdeg1);
end
sim.local.wi = wi;

% Integration weights for Lobatto integration
[xloc,wi]=lobatto_compute(sim.setup.intdeg);
sim.local.wi_quad = wi;
sim.local.xloc_quad = xloc;
sim.global.wi_quad = kron(eye(sim.setup.nx), sim.local.wi_quad);

% Polynomial values at integration points
sim.local.polyvals = zeros(sim.setup.intdeg, sim.setup.pdeg1);
for i = 1:sim.setup.pdeg1
    sim.local.polyvals(:, i) = polyval(sim.setup.poly(:,i), xloc);
end
sim.global.polyvals = kron(eye(sim.setup.nx), sim.local.polyvals);

% row and column locations for quad integral sparse matrices
sim.int.rows = reshape(kron(1:sim.setup.nx*sim.setup.pdeg1, ones(sim.setup.pdeg1,1)),1,[]);
temp = sim.setup.pdeg1*reshape(kron(1:sim.setup.nx, ones(sim.setup.pdeg1,1)),1,[]);
cols = [];
for i = 1:sim.setup.pdeg1
    cols = [temp - i+1; cols];
end
sim.int.cols = reshape(cols,1,[]);


% Global term rearranged for integration involving psi_i or psi_j
sim.global.polyvals2 = zeros(sim.setup.nx*sim.setup.intdeg, sim.setup.nx, sim.setup.pdeg1);
for i = 1:sim.setup.pdeg1
    sim.global.polyvals2(:, :, i) = kron(eye(sim.setup.nx), sim.local.polyvals(:,i));
end

%Create the permanent local matrices
poly = sim.setup.poly;
nx = sim.setup.nx;
pdeg1 = sim.setup.pdeg+1;
np = nx*pdeg1;

sim.local.P00 = MM(1,poly,poly, sim).';

sim.local.P01 = MM(1,poly,pdiff(poly,1), sim).';

sim.local.Psi = Tensor(poly,poly,poly,[0,0,1],sim);
sim.local.Psis3 = zeros(pdeg1, pdeg1*pdeg1);
sim.local.Psis1 = zeros(pdeg1, pdeg1*pdeg1);
for i = 1:pdeg1
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    sim.local.Psis1(:,index1:index2) = squeeze(sim.local.Psi(:,i,:)).';
    sim.local.Psis3(:,index1:index2) = squeeze(sim.local.Psi(:,i,:));
end

sim.local.pl = zeros(pdeg1,1);
for i = 1:pdeg1
    sim.local.pl(i) = polyval(poly(:,i),-1);
end

sim.local.pr = zeros(pdeg1,1);
for i = 1:pdeg1
    sim.local.pr(i) = polyval(poly(:,i),1);
end

[sim.local.D, sim.local.DR, sim.local.DL] = IP(1,poly,poly);

% S is L-R
[sim.local.S, sim.local.SR, sim.local.SL] = IP(1,poly,1);
sim.local.S = -sim.local.S;

% Create global matrixes for length scaling
sim.global.lt = zeros(np);
for i = 1:sim.setup.nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    sim.global.lt(index1:index2,index1:index2) =  sim.setup.lt(i) * eye(pdeg1);
end

sim.global.G = zeros(np,1);
for i = 1:sim.setup.nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    sim.global.G(index1:index2) = sim.setup.lt(i) * MM(1, sim.sbp.G(:,i), sim.setup.poly, sim);
end

sim.global.alpha = zeros(np,1);
for i = 1:sim.setup.nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    sim.global.alpha(index1:index2) =  sim.setup.lt(i) * MM(1, sim.sbp.alpha(:,i), sim.setup.poly, sim);
end

sim.global.Vn = zeros(np, np);
sim.global.Vp = zeros(np, np);
for i = 1:nx
    index1 = (i-1)*(pdeg1)+1;
    index2 = i*pdeg1;
    sim.global.Vn(index1:index2,index1:index2) = - sim.setup.lt(i) * MM(sim.sbp.En(:,i), poly, poly, sim).';
    sim.global.Vp(index1:index2,index1:index2) =  - sim.setup.lt(i) * MM(sim.sbp.Ep(:,i), poly, poly, sim).';
end

if sim.setup.tensor
    sim.global.Psi = zeros(np, np, np);
    for i = 1:nx
        index1 = (i-1)*pdeg1+1;
        index2 = i*pdeg1;
        sim.global.Psi(index1:index2,index1:index2,index1:index2) = sim.local.Psi;
    end
else
    
    sim.global.Psis1 = sparse(np, np*pdeg1);
    sim.global.Psis3 = sparse(np, np*pdeg1);
    
    for i = 1:nx
        index1 = (i-1)*pdeg1+1;
        index2 = i*pdeg1;
        
        index3 = (i-1)*pdeg1^2+1;
        index4 = i*pdeg1^2;
        
        sim.global.Psis1(index1:index2,index3:index4) = sim.local.Psis1;
        sim.global.Psis3(index1:index2,index3:index4) = sim.local.Psis3;
    end
end



sim.global.Psialphas1 = sparse(np, np*pdeg1);
sim.global.Psialphas3 = sparse(np, np*pdeg1);
for i = 1:nx
    alpha_ini2_poly = 0*poly;
    for j = 1:pdeg1
        alpha_ini2_poly(:, j) = polymatch(polytimes(poly(:,j),polytimes(sim.sbp.alpha(:,i), sim.sbp.ini2(:,i))), pdeg1);
    end
    
    Psialpha = Tensor(alpha_ini2_poly,poly,poly,[0,0,0],sim);
    sim.local.Psialphas3 = zeros(pdeg1, pdeg1*pdeg1);
    sim.local.Psialphas1 = zeros(pdeg1, pdeg1*pdeg1);
    for j = 1:pdeg1
        index1 = (j-1)*pdeg1+1;
        index2 = j*pdeg1;
        sim.local.Psialphas1(:,index1:index2) = squeeze(Psialpha(:,j,:)).';
        sim.local.Psialphas3(:,index1:index2) = squeeze(Psialpha(:,j,:));
    end
    
    
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    
    index3 = (i-1)*pdeg1^2+1;
    index4 = i*pdeg1^2;
    
    sim.global.Psialphas1(index1:index2,index3:index4) = sim.setup.lt(i) * sim.local.Psialphas1;
    sim.global.Psialphas3(index1:index2,index3:index4) = sim.setup.lt(i) * sim.local.Psialphas3;
end


% Create global versions
sim.global.P00 =  sim.global.lt * kron(eye(nx), sim.local.P00);

sim.global.P01 = kron(eye(nx), sim.local.P01);

sim.global.DR = kron(eye(nx), sim.local.DR.');
sim.global.DL = kron(eye(nx), sim.local.DL.');
sim.global.D1 = sim.global.DR-sim.global.DL;

sim.global.SR = padarray(kron(eye(nx-1), sim.local.SR), [0,pdeg1],  0, 'post');
sim.global.SL = padarray(kron(eye(nx-1), sim.local.SL), [0,pdeg1],  0, 'pre' );

sim.global.S_psi_n = sim.global.DeltaChi * sim.global.SL - sim.global.SR;
sim.global.S_psi_p = sim.global.DeltaChiEg * sim.global.SL - sim.global.SR;
sim.global.S_psi = sim.global.SL - sim.global.SR;

sim.global.Pn = zeros(np, np);
for i = 1:nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    sim.global.Pn(index1:index2,index1:index2) =  sim.setup.lt(i) *  MM(sim.sbp.c_n(:,i), poly, poly, sim).';
end

sim.global.Pp = zeros(np, np);
for i = 1:nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    sim.global.Pp(index1:index2,index1:index2) =  sim.setup.lt(i) *  MM(sim.sbp.c_p(:,i), poly, poly, sim).';
end

sim.global.Pphi = zeros(np, np);
for i = 1:nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    sim.global.Pphi(index1:index2,index1:index2) =  sim.setup.lt(i) * MM(sim.sbp.c_phi(:,i), poly, poly, sim).';
end


end




