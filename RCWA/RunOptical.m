

function sim = RunOptical(sim)


testval = [inf, inf];
sim.results.JscOpt = cell(sim.setup.nlambda,1);

for ll = 1 : sim.setup.nlambda
    
    % Reset parameters for Nt optimisation
    sim.results.success = false;
    sim.setup.newNt = sim.setup.Nt;
    if sim.setup.CheckNtConvergence
        sim.results.Ntprev = nan;
        sim.results.Nttest(ll) = nan;
        sim.results.failcount(ll) = 0;
        sim.results.succount(ll) = 0;
        sim.results.failcount(ll) = 0;
        sim.results.predictedNtprev(ll) = nan;
        sim.results.predictedNt(ll) = nan;
        
    end
    
    while ~sim.results.success
        % We now want epsf(z, l)
        sim.setup.Nt = sim.setup.newNt;
        sim = BuildFourier(sim);
        
        % Set permittivity for current wavelength
        eps_struct1 = sim.zf.eps_struct1{sim.setup.Nt+1}; % Primary permittivity (i.e. metal)
        eps_struct2 = sim.zf.eps_struct2{sim.setup.Nt+1}; % Secondarty grating permittivity
        eps1 = sim.zl.eps(:,ll);
        eps2 = sim.zl.eps2(:,ll);
        sim.zfl.eps{ll} = diag(eps1) * eps_struct1+diag(eps2) * eps_struct2;
        
        % Set inverse permittivity for current wavelength
        ieps1 = 1./sim.zl.eps(:,ll);
        ieps1(ieps1 == Inf) = 0;
        ieps2 = 1./sim.zl.eps2(:,ll);
        ieps2(ieps2 == Inf) = 0;
        sim.zfl.ieps{ll} = diag(ieps1) * eps_struct1 + diag(ieps2) * eps_struct2;
        
        
        % Create zxl permittivity
        % Set to eps2, then correct using grating shape
        Nx = sim.setup.Nx;
        sim.zx.eps = kron(eps2,ones(1, Nx));
        
        % Find upper and lower limit of x for central metal (eps1) section
        locations = floor(sim.z.zetalist*(Nx-1)+1+Nx/2);
        for i = 1:floor(size(locations,1)/2)
            for j = 1:size(locations,2)
                xdielower = locations(2*i-1,j);
                xdieupper = locations(2*i,j);
                sim.zx.eps(j, xdielower:xdieupper) = eps1(j);
            end
            
            % sim.zxl.epsif = InverseFourier(sim.setup.nmx', sim.zfl.eps{ll}, sim.input.nmLx, 0);
        end
        
        % Variable theta not implemented
        for tt = 1:1
            if sim.setup.pol == 0
                for i = 1:2
                    sim.setup.pol = i;
                    sim = RCWA(ll, tt, sim);
                end
                sim.setup.pol = 0;
            else
                sim = RCWA(ll, tt, sim);
            end
        end
        
        
        % Calculate absorption
        sim  = AnalyseAbs(sim, ll);
        
        % Calculate convergence indicator
        testval = fliplr(testval);
        testval(2) = sim.results.mAicm2inmJOpt(ll)*10^9;
        
        if sim.setup.CheckNtConvergence
            sim = checkConverge(sim, testval, ll);
        else
            sim.results.success = 1;
        end
    end
    
    % Store convergence indicator
    sim.results.testval(ll) = testval(2);
    testval = [inf, inf];
    sim.results.Nt(ll) = sim.setup.Nt;

    if sim.setup.CheckNtConvergence
        sim.setup.Nt = max(sim.setup.Nt-8, sim.setup.minNt);
    end
    
end

% Calculate absoption, G etc.
sim  = AnalyseAbs(sim);
if sim.setup.optplot == 1
    figure()
    for ll = 1:sim.setup.nlambda
        plot(sim.results.JscOpt{ll})
        hold on
    end
end

if sim.setup.pol == 0 && sim.setup.RebuildFields
    sim.zx.Ex = 1/sqrt(2) * sim.zx.Ex;
    sim.zx.Ey = 1/sqrt(2) * sim.zx.Ey;
    sim.zx.Ez = 1/sqrt(2) * sim.zx.Ez;
    
    %sim.zx.Hx = 1/sqrt(2) * sim.zx.Hx;
    sim.zx.Hy = 1/sqrt(2) * sim.zx.Hy;
    %sim.zx.Hz = 1/sqrt(2) * sim.zx.Hz;
end

if sim.setup.optplot == 1
    if length(sim.setup.nmlambda) > 1
        figure(1);
        hold off
        plot(sim.setup.nmlambda,sim.results.sA, 'b:', sim.setup.nmlambda, sim.results.pA, 'r:', sim.setup.nmlambda,sim.results.Ajunc_s, 'b', sim.setup.nmlambda, sim.results.Ajunc_p, 'r', sim.setup.nmlambda, sim.results.Ajunc, 'g', sim.setup.nmlambda, sim.results.Arest, 'r:');
        hold on
        plot(sim.setup.nmlambda, sim.results.Nt/max(sim.results.Nt))
        legend('Total s', 'Total p', 'Junc s', 'Junc p', 'Junction unpol');
    end
    
    if (sim.setup.pol == 0 || sim.setup.pol == 2)&& sim.setup.Eplot > 0
        figure();
        colormap(jet);
        %     % Different graphing options
        contourf(sim.setup.nmx, -sim.setup.nmz, abs(sim.zx.Ey), 200,'linestyle','none');
        title('Field: Ey');
        caxis([0,inf])
        axis equal;
        xlim([sim.setup.nmx(1), sim.setup.nmx(end)]);
    end
    
    if (sim.setup.pol == 0 || sim.setup.pol == 1) && sim.setup.Hplot>0
        figure();
        colormap(jet);
        %     % Different graphing options
        contourf(sim.setup.nmx, -sim.setup.nmz, abs(sim.zx.Hy), 200,'linestyle','none');
        title('Field: Hy');
        caxis([0,inf])
        axis equal;
        xlim([sim.setup.nmx(1), sim.setup.nmx(end)]);
    end
    
    if sim.setup.epsplot>0
        figure();
        colormap(jet);
        %     % Different graphing options
        contourf(sim.setup.nmx, -sim.setup.nmz, abs(sim.zx.eps), 200,'linestyle','none');
        title('Field: abs(eps)');
        caxis([0,inf])
        axis equal;
        xlim([sim.setup.nmx(1), sim.setup.nmx(end)]);
    end
end

end
