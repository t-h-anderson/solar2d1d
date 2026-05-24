
function sim = BuildJV(sim)
% Buld JV curve, running in both positive and negative steps. Take minimum
% of each direction to produce most converged values.

% Initialise
sim.results.JVforward = [];
sim.results.PVforward = [];
J=0;
sim.setup.Vext = 0;
Jsign = 1; % Toggle for sign change in J
P=0;
ptoggle = 0; % Save sim near maximum P

if ~isfield(sim.input,'mAicm2JOpt')
    sim.input.mAicm2JOpt = inf;
end

% Build upper quadrant of JV curve
while  J >= 0 && J <=sim.input.mAicm2JOpt
    
    % Converge
    sim = Converge(sim);
    
    if sim.results.fail 
        if sim.setup.dV > sim.setup.mindV
            dV = sim.setup.dV;
            sim = simsave;
            sim.setup.dV = dV/2;
        else
            break;
        end
    else
        
        % Calculate J from coeff list
        Jspace = coeffs2space(sim.coeffs.Jn + sim.coeffs.Jp, sim);
        
        % Take J to be average J in device
        J = Jsign * sim.phys.mAicm2Js * sum(Jspace(:,2))/size(Jspace,1);
        
        % If first J, check sign and reverse if necessary to make positive
        if isempty(sim.results.JVforward) && J < 0
            Jsign = -1;
            J = Jsign * J;
        %elseif (size(sim.results.JVforward, 2) ==1 && J < 0) 
        % For Faiz, untested, may help if inbuilt bias is tiny.
        %    dV = sim.setup.dV;
        %    sim = simsave;
        %    sim.setup.dV = dV/2;
        %    break;
        end
        
        % Display results
        if sim.setup.SolverSteps
            fprintf('Vext: %.2f, J: %.2f, Loop: %d, Abs: %.1g, Rel: %.2g \n', sim.phys.VVs * sim.setup.Vext, J, sim.setup.counter, sim.setup.abschange, sim.setup.relchange);
        end
        
        % Append results
        sim.results.JVforward = [sim.results.JVforward; [- sim.phys.VVs * sim.setup.Vext, J]];
        sim.results.PVforward = [sim.results.PVforward; [- sim.phys.VVs * sim.setup.Vext, - sim.phys.VVs * sim.setup.Vext*J]];
        
        
        % If first local maximum power, save point
        if  abs(sim.setup.Vext* J) >= P && ptoggle == 0
            simsave = sim;
            simsave.results.P = abs(sim.setup.Vext* J) * sim.phys.VVs;
            simsave.results.Vmax0 = abs(sim.phys.VVs * sim.setup.Vext);
        else
            ptoggle= 1;
        end
        
        % Increase loop variables
        P = max(P, abs(sim.setup.Vext* J));
    end
    sim.setup.Vext = sim.setup.Vext + Jsign * sim.setup.dV;
end

if sim.results.fail
    JV = [0,0];
else
    
    % Rerun calculation backward
    if sim.setup.rerunbackward
        P2=0;
        sim.results.JVbackward = [];
        for i = size(sim.results.JVforward, 1):-1:1
            
            sim.setup.Vext = -sim.results.JVforward(i,1)/sim.phys.VVs;
            
            % Converge
            sim = Converge(sim);
            
            % Calculate J from coeff list
            Jspace = coeffs2space(sim.coeffs.Jn + sim.coeffs.Jp, sim);
            
            % Take J to be average J in device
            J = Jsign * sim.phys.mAicm2Js * sum(Jspace(:,2))/size(Jspace,1);
            
            % Display results
            if sim.setup.SolverSteps
                fprintf('Vext: %.2f, J: %.2f, Loop: %d, Abs: %.1g, Rel: %.2g \n', sim.phys.VVs * sim.setup.Vext, J, sim.setup.counter, sim.setup.abschange, sim.setup.relchange);
            end
            
            % Append results
            sim.results.JVbackward = [sim.results.JVbackward; [- sim.phys.VVs * sim.setup.Vext, J]];
            
            % If first local maximum power, save point
            if  abs(sim.setup.Vext* J) >= P2 && ptoggle == 0
                simsave = sim;
            else
                ptoggle= 1;
            end
            
            % Increase loop variables
            P2 = max(P2, abs(sim.setup.Vext* J));
        end
        
        sim.results.JV = [sim.results.JVforward(:,1), min(sim.results.JVforward(:,2), flipud(sim.results.JVbackward(:,2)))];
        
        if sim.setup.SolverSteps
            Jdiff = 2*max(abs(sim.results.JVforward(:,2) - flipud(sim.results.JVbackward(:,2)))./(sim.results.JVforward(:,2) + flipud(sim.results.JVbackward(:,2))));
            if Jdiff > 0.01
                disp('Warning, JV hysterisis.')
            end
            if sim.setup.plotsolve>0
                plot(-sim.results.JVforward(:,1), sim.results.JVforward(:,2),-sim.results.JVforward(:,1),flipud(sim.results.JVbackward(:,2)))
            end
            
        end
    else
        sim.results.JV = sim.results.JVforward;
    end
    % resolve maximum JVP
    JV = sim.results.JV;
    P = abs(JV(:,1)).*JV(:,2);
    
    % Find maximum P
    [Vsort, ind] = max(P);
    
    Ppoints= size(JV, 1);
    
    % Create 5 point stencil about maximum power point
    if Ppoints > 1
        if ind > 1 && ind+1 <= Ppoints % Best case
            JVpeak = [JV(ind-1,1), 0.5*(JV(ind-1,1)+JV(ind,1)), JV(ind,1), 0.5*(JV(ind,1) + JV(ind+1,1)), JV(ind+1,1); ...
                JV(ind-1,2), 0, JV(ind,2), 0, JV(ind+1,2)];
        elseif ind ==1 % ind = 1
            JVpeak = [JV(1,1), 0.75*JV(1,1)+0.25*JV(2,1), 0.5*JV(1,1)+0.5*JV(2,1), 0.25*JV(1,1)+0.75*JV(2,1), JV(2,1); ...
                JV(1,2), 0, 0, 0, JV(2,2)];
        else
            JVpeak = [JV(end-1,1), 0.75*JV(end-1,1)+0.25*JV(end,1), 0.5*JV(end-1,1)+0.5*JV(end,1), 0.25*JV(end-1,1)+0.75*JV(end,1), JV(end,1); ...
                JV(end-1,2), 0, 0, 0, JV(end,2)];
        end
    else
        
        JV = [0,0];
    end
    
    
    % Reload saved simulation point
    temp = sim.results;
    sim = simsave;
    sim.results = temp;
    i = 0;
    Jrel = inf;
    
    % Try to improve the MPP 10 times
    while i < sim.setup.PVrefinemax && Jrel > sim.setup.PVtol
        
        i = i + 1;
        if J >=sim.input.mAicm2JOpt
            break;
        end
        
        % Fill in the blanks
        Vtry = JVpeak(1,JVpeak(2,:)==0);
        for j = 1:length(Vtry)
            sim.setup.Vext = -Vtry(j)/sim.phys.VVs;
            
            % Converge
            newsim = Converge(sim);
            
            
            % Calculate J from coeff list
            Jspace = coeffs2space(newsim.coeffs.Jn + newsim.coeffs.Jp, newsim);
            
            % Take J to be average J in cell
            J = newsim.phys.mAicm2Js * Jsign * sum(Jspace(:,2))/size(Jspace,1);
            
            % Display results
            if newsim.setup.SolverSteps
                fprintf('Vext: %.4f, J: %.2f, Loop: %d, Abs: %.1g, Rel: %.2g \n', newsim.phys.VVs * newsim.setup.Vext, J, newsim.setup.counter, newsim.setup.abschange, newsim.setup.relchange);
            end
            
            sim.results.JV = [sim.results.JV; [- newsim.phys.VVs * newsim.setup.Vext, J]];
            JVpeak(2, JVpeak(1,:) == Vtry(j)) = J;
        end
        
        
        % Check to make sure solution is convex
        %maxmin = sign(conv(conv(abs(JVpeak(1,:).*JVpeak(2,:)),[1,-1],'valid'),[-1,1],'valid'));
        %if abs(sum(maxmin)) < 2
        %    disp('JV curve failed to resolve');
        %    break
        %end
        
        Jrel = (max(JVpeak(2,:)) - min(JVpeak(2,:)))/ min(JVpeak(2,:));
        P = abs(JVpeak(1,:).*JVpeak(2,:));
        [Psort,ind] = sort(P, 'descend');
        Ppoints= size(JVpeak, 2);
        
        % Create 5 point stencil about maximum power point
        ind = ind(1);
        if Ppoints > 1
            if ind > 1 && ind+1 <= Ppoints % Best case
                JVpeak = [JVpeak(1, ind-1), 0.5*(JVpeak(1, ind-1)+JVpeak(1, ind,1)), JVpeak(1, ind), 0.5*(JVpeak(1, ind) + JVpeak(1, ind+1)), JVpeak(1, ind+1); ...
                    JVpeak(2, ind-1), 0, JVpeak(2, ind), 0, JVpeak(2, ind+1)];
            elseif ind ==1 % ind = 1
                JVpeak = [JVpeak(1,1), 0.75*JVpeak(1,1)+0.25*JVpeak(1, 2), 0.5*JVpeak(1,1)+0.5*JVpeak(1, 2), 0.25*JVpeak(1,1)+0.75*JVpeak(1, 2), JVpeak(1, 2); ...
                    JVpeak(2, 1), 0, 0, 0, JVpeak(2,2)];
            else
                JVpeak = [JVpeak(1, end-1), 0.75*JVpeak(1, end-1)+0.25*JVpeak(1, end), 0.5*JVpeak(1, end-1)+0.5*JVpeak(1, end), 0.25*JVpeak(1, end-1)+0.75*JVpeak(1, end), JVpeak(1, end); ...
                    JVpeak(2, end-1), 0, 0, 0, JVpeak(2, end)];
            end
        else
            error('Only one power point computed')
        end
        
    end
end
sim.results.J = sim.results.JV(:,2);
sim.results.V = sim.results.JV(:,1);
sim.results.P = sim.results.J.*abs(sim.results.V);

[sim.results.V,ind] = sort(abs(sim.results.V));
sim.results.J = sim.results.J(ind);
sim.results.P = sim.results.P(ind);
end
