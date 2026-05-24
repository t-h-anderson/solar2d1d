
function result = intdnmlambda(var, sim)
% Integrate over wavelength
if length(size(var)) == 3
    if sim.setup.nlambda > 1
        result = (sum(var, 3) - 0.5*var(:,:,1) - 0.5*var(:,:,end)) * sim.setup.nmdlambda;
    else
        result = var(:,:,1);
    end
elseif length(size(var)) == 2
    if sim.setup.nlambda > 1
        result = (sum(var, 2) + 0.5*var(:,1) + 0.5*var(:,end)) * sim.setup.nmdlambda;
    else
        result = var(:,1);
    end
else
    error('intdnmlambda not defined for lengths other than 2 and 3');
end


end

