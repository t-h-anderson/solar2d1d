

function sim = BuildFourier(sim)

Nt = sim.setup.Nt;

if ~isfield(sim,'zf')
    % Fourier series for eps

    sim.zf.eps_struct1{Nt+1} = zeros(sim.setup.Nz, 4*Nt+1);
    sim.zf.eps_struct2{Nt+1} = zeros(sim.setup.Nz, 4*Nt+1);
end

% Load previous fourier, or create new
Ntmax = length(sim.zf.eps_struct1)-1;
if Nt <= Ntmax
    epsf1 = sim.zf.eps_struct1{Nt+1};
    epsf2 = sim.zf.eps_struct2{Nt+1};
end

if ~exist('epsf1', 'var')
    epsf1 = zeros(sim.setup.Nz, 4*Nt+1);
    epsf2 = zeros(sim.setup.Nz, 4*Nt+1);
end

if isempty(epsf1)
    epsf1 = zeros(sim.setup.Nz, 4*Nt+1);
    epsf2 = zeros(sim.setup.Nz, 4*Nt+1);
end

% If new, calculate epsf1 and epsf2
if epsf1(1,2*Nt+1) ~= 1
    
    epsf1(:, 2*Nt +1) = 1;
    
    ngsec = size(sim.z.zetalist,1)-1;
    for i = 1:sim.setup.Nz
        for n = -2*sim.setup.Nt:2*sim.setup.Nt     
            if n ~= 0
                fstructure1 = 0;
                for j = 1:ceil(ngsec/2)
                    fstructure1 = fstructure1 + exp(2*n.*pi.*1.i.*sim.z.zetalist(2*j,i)) - exp(2*n.*pi.*1.i.*sim.z.zetalist(2*j-1,i));
                end
                fstructure1 = fstructure1/(2.i*n*pi);
                
                fstructure2 = 0;
                for j = 1:floor(ngsec/2)
                    fstructure2 = fstructure2 + exp(2*n.*pi.*1.i.*sim.z.zetalist(2*j+1,i)) - exp(2*n.*pi.*1.i.*sim.z.zetalist(2*j,i));
                end
                fstructure2 = fstructure2/(2.i*n*pi);
            else
                
                fstructure1 = 0;
                for j = 1:ceil(ngsec/2)
                    fstructure1 = fstructure1 + sim.z.zetalist(2*j,i) - sim.z.zetalist(2*j-1,i);
                end
                
                fstructure2 = 0;
                for j = 1:floor(ngsec/2)
                    fstructure2 = fstructure2 + sim.z.zetalist(2*j+1,i) - sim.z.zetalist(2*j,i);
                end
                
            end
            epsf1(i, n + 2*Nt+1) = fstructure1;
            epsf2(i, n + 2*Nt+1) = fstructure2;
        end
    end
end

sim.zf.eps_struct1{Nt+1} = epsf1;
sim.zf.eps_struct2{Nt+1} = epsf2;
end





