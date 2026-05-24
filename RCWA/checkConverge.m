function sim = checkConverge(sim, fom, ll)

% Save old values
sim.results.Nttestprev(ll) = sim.results.Nttest(ll);

if isinf(sim.results.Nttestprev(ll))
    sim.results.Nttestprev(ll) = nan;
end

%Allow plot of how Jsc changes with Nt
sim.results.JscOpt{ll} = [sim.results.JscOptll{ll}, fom(2)];

% Test change in fom
if sum(fom)~=0
    Nttest = 2*sim.setup.nmdlambda * abs(fom(2) - fom(1));%/sum(fom);
    sim.results.Nttest(ll) = Nttest;
else
    sim.results.Nttest(ll) = 0;
end

predictedNt = sim.setup.Nt+2;

% If test if greater than tolerance
if  sim.results.Nttest(ll) > sim.setup.Nttol
    predtest = nan;
    % If test has decreased
    if sim.results.Nttest(ll) < sim.results.Nttestprev(ll)
        % Predict required Nt
        % Fit two points to exponential
        sim.results.predictedNtprev(ll) = sim.results.predictedNt(ll);

        fitfn = fit([sim.setup.Nt; sim.results.Ntprev], [sim.results.Nttest(ll); sim.results.Nttestprev(ll)], 'exp1');
        sim.results.predictedNt(ll) = ceil(log(sim.setup.Nttol/fitfn.a)/fitfn.b);
                       
        % If subsequent predictions are similar, take the larger
        predtest = abs(sim.results.predictedNt(ll)-sim.results.predictedNtprev(ll))/mean([sim.results.predictedNt(ll),sim.results.predictedNtprev(ll)]);
        if predtest < sim.setup.ptol
            testNt = 2*ceil(mean([sim.results.predictedNt(ll),sim.results.predictedNtprev(ll)]-1)/2)+1;
            predictedNt = min(max(testNt,predictedNt),2*predictedNt);
        end
    else
        % Else mark bad point
        sim.results.failcount(ll) = sim.results.failcount(ll) + 1;
    end
    
%     if loc.plotting == 1
         fprintf('No! \tlambda: %.3g,\tNt: %g,\tfom:%.3g,\ttest:%.3g,\tpredtest:%.2g,\tPredicted Nt:%g,\tPredicted Nt old:%g\n', ...
             sim.setup.nmlambda(ll), sim.setup.Nt, fom(2), sim.results.Nttest(ll), predtest, sim.results.predictedNt(ll), sim.results.predictedNtprev(ll));
         sim.succount = 0;
%     end
% elseif loc.Nttest == inf
%     sim.succount = sim.nosuccreq;
elseif sim.results.Nttest(ll) < sim.setup.Nttol_lower
%     if loc.plotting == 1
         fprintf('WOW!\tlambda: %.3g,\tNt: %g,\tfom:%.3f,\ttest:%.3g\n', sim.setup.nmlambda(ll), sim.setup.Nt, fom(2), sim.results.Nttest(ll));
%     end
     sim.results.succount(ll) = sim.setup.nosuccreq;
else
%     if loc.plotting == 1
         fprintf('OK! \tlambda: %.3g,\tNt: %g,\tfom:%.3g,\ttest:%.3g\n', sim.setup.nmlambda(ll), sim.setup.Nt, fom(2), sim.results.Nttest(ll));
%     end
    sim.results.succount(ll) = sim.results.succount(ll) + 1;
end

if sim.results.succount(ll) == sim.setup.nosuccreq
    sim.results.success = true;
end

if sim.results.failcount(ll) > sim.setup.faillimit
%     if loc.plotting == 1
         fprintf('Warning! Max convergence failed %.0f times\n', sim.setup.faillimit);
%     end
    sim.success = true;
end

if sim.setup.Nt >= sim.setup.maxNt
    if sim.setup.optplot == 1
        fprintf('Warning! Max Nt = %.0f reached\n',sim.setup.maxNt);
    end
    sim.success = true;
end

% Update predicted Nt
sim.results.Ntprev = sim.setup.Nt;
sim.setup.newNt = min(max(predictedNt-2, sim.setup.Nt)+2, sim.setup.maxNt);

end
