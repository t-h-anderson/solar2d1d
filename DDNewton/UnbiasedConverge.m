
function sim = UnbiasedConverge(sim)

loopcount = 1; % Counter for max loops allowed
negfail = 0;
negcount = 0;
maxnegloops = 10; % Max subs no. steps where there are negative carriers

while loopcount <= sim.setup.maxloop && sim.results.fail == 0
        
    % Increase loopcount
    loopcount = loopcount + 1;
    
    % Initialise bc and hybridization
    sim = set_t(sim);
    sim = setboundaries(sim);

    % Perform Newton Iteration
    [sim, abschange, relchange] = NewtonSolve(sim);
    sim.setup.time1 = toc(sim.setup.start_time); %replaced start_electrical to start_time:Faiz
    
    % Display results
    if sim.setup.SolverSteps
        fprintf('Damping: %d, Iteration: %d, Abs: %.1g, Rel: %.2g \n', sim.setup.damping, sim.setup.counter, abschange, relchange);
    end
    
    % Kill if any NaNs
    if max(isnan(sim.coeffs.n_hat + sim.coeffs.p_hat)) >0 || max(sim.coeffs.n_hat + sim.coeffs.p_hat) > 10^10
        disp('NaN found, increasing initial damping')
        sim.results.fail = 1;
    end
    
    % Check tolerances
    if (relchange < sim.setup.reltol && abschange < sim.setup.abstol)
        if sim.setup.damping <1
            % Reduce damping
            loopcount = 1; % Reset counter
            
            if sim.setup.time1 < sim.setup.maxtime
                sim.setup.damping = min(sim.setup.damprelax * sim.setup.damping, 1); % Reduce damping
            else
                sim.setup.damping = 1;
            end
            sim = initializeDevice(sim);
            sim = BCond(sim);
        else
            if min([sim.coeffs.n_hat; sim.coeffs.p_hat]) > 0 || negfail > 3
                % Success!
                break;
            else
                % Negative carriers, try one more loop
                negfail = negfail + 1;
            end
        end
    end
    
    % Draw graphs
    if sim.setup.plotsolve > 0 && mod(loopcount, sim.setup.plotsolve) == 0
        displayResults(sim);
    end
    
    
    % If loopcount has maxed out, damping not finished and carriers are positive, relax and warn
    if loopcount >= sim.setup.maxloop && sim.setup.damping < 1 && min([sim.coeffs.n_hat; sim.coeffs.p_hat]) >= 0
        if sim.setup.SolverSteps
            disp('Warning, simulation not converged within maxloop steps. Damping relaxed anyway. Good luck!');
        end
        loopcount = 1; % Reset counter
        sim.setup.damping = min(sim.setup.damprelax * sim.setup.damping, 1); % Reduce damping
        sim = initializeDevice(sim);
        sim = BCond(sim);
    end
    
    % Check negative iterations
    if min([sim.coeffs.n_hat; sim.coeffs.p_hat]) < 0
        negcount = negcount + 1;
    else
        negcount = 0;
    end
    
    if negcount > maxnegloops
        sim.results.fail = 1;
    end
    
end

% If skeleton carriers are negative, fail point
if min([sim.coeffs.n_hat; sim.coeffs.p_hat]) < 0
    if sim.setup.plotsolve > 0
        displayResults(sim);
    end
    
    % Flag where negative
    if sim.setup.refine_on_negative
        npneg = sign((sim.coeffs.n_hat < 0) + (sim.coeffs.p_hat < 0));
        sim.setup.refinesecs = sign(conv([1,1],npneg>sim.setup.Jtol));
    end
    
    if sim.setup.warn
        disp('Warning, negative carriers.');
    end
    
    sim.results.fail = 1;
end

if ~sim.results.fail
    
    % Test how smooth the current skeleton is
    J = coeffs2sbp(sim.coeffs.Jn+sim.coeffs.Jp,sim);
    Jconst = abs(J(end));
    Jnonconst = sum(abs(J(1:end-1,:)),1);
    Jnoise = abs(Jnonconst./Jconst);
    [~, worst]=sort(Jnoise,'descend');
    Jnoise(worst(ceil(sim.setup.nx*sim.setup.maxJnoiserefinefrac):end)) = 0;
    
    % If noise is bad, send refine locations
    Jtest = max(Jnoise) ;
    if  Jtest > sim.setup.Jtol
        sim.results.fail = 1;
        
        if sim.setup.plotsolve > 0
            displayResults(sim);
        end
        
        if sim.setup.warn
            disp('noisy J found')
        end
        
        % Flags for mesh refinement - needs some work to aid convergence
        temp=sign(conv([1,1,1], Jnoise > sim.setup.Jtol));
        sim.setup.refinesecs = temp(1:end-1);
    end  
end




end

