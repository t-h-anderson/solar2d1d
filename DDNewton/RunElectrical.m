

function sim = RunElectrical(sim)

trydamping = sim.setup.trydamping;

sim = initializeDevice(sim);
sim = initializeSim(sim);

damp0 = sim.setup.dampmin;

sim.setup.time1 = 0;
while trydamping >= 0 && sim.setup.time1 < sim.setup.maxtime
    
    % set boundaries conditions
    sim = BCond(sim);
    sim = SetInitialConditions(sim);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Final Setup        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sim.setup.plotsolve > 0
        displayResults(sim);
    end
    
    sim.results.fail = 0;
    sim.setup.refinesecs = 0;
    
    sim = UnbiasedConverge(sim);
    
    if sim.setup.plotsolve > 0
        displayResults(sim);
    end
    
    % Break if succeeded
    if ~sim.results.fail
        break;
    end
    
    trydamping = trydamping - 1;
    
    if trydamping>=0
        % Tweak parameters
        if sim.setup.pdeg < sim.setup.maxpdeg
            sim.setup.pdeg = sim.setup.pdeg + 1;
        else
            if sum(sim.setup.refinesecs) > 0
                %sim.setup.refinesecs = 0*sim.setup.mesh + 1;
            end
        end
        
        damp0=damp0/sim.setup.dampingincrease;
        sim.setup.dampmin = damp0;
        
        % Reset Sim - apply refinesecs
        sim = initializeDevice(sim);
        sim = initializeSim(sim);
        
        if sim.setup.plotsolve > 0
            displayResults(sim);
        end
    end
    
    % Update time to terminate after limit reached
    sim.results.time1 = toc(sim.setup.start_time);%replaced start_electrical to start_time:Faiz
    
end

if sim.setup.ExportSCdetails
    exportResults(sim);
end

if ~sim.results.fail
    % Ramp J here
    sim = BuildJV(sim);
    
    if sim.setup.plotsolve > 0
        displayResults(sim);
    end
    
    % Calculate the short circuit current
    sim.results.mAicm2Jsc = sim.results.JV(1,2);
    
    if length(sim.input.icm3isG) == 1
        sim.results.JscOpt = 1000 * cm_from_nm(sim.setup.nmLz) * sim.phys.Cq * sim.input.icm3isG;
    else
        GFn = str2func(sim.input.icm3isG);
        z = 0:0.01:1;
        sim.results.JscOpt = 1000 * cm_from_nm(sim.setup.nmLz) * sim.phys.Cq * sum(GFn(z))/length(z);
    end
    
    % Calculate the open circuit voltage
    JV = sim.results.JV;
    mAicm2J = JV(:,2);
    VV = abs(JV(:,1));
    [VV, order] = sort(VV);
    mAicm2J = mAicm2J(order);
    dJ = conv(mAicm2J,[1,-1],'valid');
    dV = conv(VV,[1,-1],'valid');
    m = dV./dJ;
    inter = [0; VV(1:end-1) - mAicm2J(1:end-1).* m];
    inter = inter.*~isinf(inter);
    sim.results.VVOC = abs(min(inter(VV>inter)));
    
    
    % Calculate the maximum power point
    Wim2P = 10*VV.*mAicm2J;
    [Wim2Pmax, Pmaxloc] = max(Wim2P);
    sim.results.Wim2Pmax = Wim2Pmax;
    sim.results.VVmax = VV(Pmaxloc);
    sim.results.mAcm2Jmax = mAicm2J(Pmaxloc);
    sim.results.FF = Wim2Pmax/(10*sim.results.VVOC*sim.results.mAicm2Jsc);
    
else
    sim.results.JV = [nan,nan];
    sim.results.VVOC = nan;
    sim.results.Wim2Pmax = nan;
    sim.results.VVmax = nan;
    sim.results.mAcm2Jmax = nan;
    sim.results.FF = nan;
end

sim.results.time1 = toc(sim.setup.start_time);%replaced start_electrical to start_time:Faiz

end
