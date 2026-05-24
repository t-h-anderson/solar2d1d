
function sim = Converge(sim)
% Iterate till solution has converged
loopcount = 1;
while loopcount <= sim.setup.maxloop
    
    % Set up simulation
    sim = BCond(sim);
    sim = set_t(sim);
    sim = setboundaries(sim);
    
    % Iterate
    
    [newsim, abschange, relchange] = NewtonSolve(sim);
    
    % Test how smooth the current skeleton is
    J = coeffs2sbp(newsim.coeffs.Jn+newsim.coeffs.Jp, newsim);
    Jconst = abs(J(end));
    Jnonconst = sum(abs(J(1:end-1,:)),1);
    Jnoise = abs(Jnonconst./Jconst);
    [~, worst]=sort(Jnoise,'descend');
    Jnoise(worst(ceil(newsim.setup.nx*newsim.setup.maxJnoiserefinefrac):end)) = 0;
    
    % If noise is bad, send refine locations
    npneg = sum(sign((sim.coeffs.n_hat < 0) + (sim.coeffs.p_hat < 0)));
    Jtest = max(Jnoise) ;
    if  npneg>0 %|| Jtest > newsim.setup.Jtol
        sim.results.fail = 1;
        break;
    else
        sim = newsim;
        
        sim.setup.abschange = abschange;
        sim.setup.relchange = relchange;
        % Break if acceptable
        if (relchange < sim.setup.reltol && abschange < sim.setup.abstol)
            break;
        end
        
        % Display results
        if sim.setup.plotsolve > 0 && mod(loopcount, sim.setup.plotsolve) == 0
            displayResults(sim);
            drawnow;
        end
    end
    loopcount = loopcount + 1;
end
end

