function [poly] = space2sbp(input, sim)
% Takes an input described sectionwise (e.g. for the i-layer) and turns it
% into a standard piecewise polynomial on the mesh.

setup = sim.setup;
phys = sim.phys;

% Allows choice of deg of poly
if ~isempty(varargin)
    pdeg1 = varargin{1};
else
    pdeg1 = setup.pdeg1;
end

% If values then 
if isa(input,'double')
    poly = zeros(pdeg1, setup.nx);
    
    sec1 = round((setup.nx* phys.Sec/phys.Lx)*tril(ones(phys.nsec))');
    sec0 = ([0, sec1(1:end-1)])+1;
    
    for i = 1:length(input)
        poly(end, sec0(i):sec1(i)) = input(i);
    end
    
elseif isa(input,'function_handle')
    poly = func2poly(input, setup.pdeg1, setup.mesh);
else
    class(input)
    
    error('Input must be a list of doubles, or an anonymous function string');
end

end

