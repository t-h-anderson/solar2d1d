function [poly] = input2sbp(input, field, scale, sim)
% Takes an input described globally and turns it
% into a standard piecewise polynomial on the mesh.

% Load the input
input = input.(field);

% If characters, split at '&' and make cells
if isa(input, 'char')
    input = strsplit(input,'&');
    for i=1:length(input)
        check = input{i};
        if ~strcmp(check(1),'@')
            input{i} = str2double(check);
        end
    end
end

% Check length of input
linput  = length(input);
if linput ~= 1 && linput ~= sim.material.nsec
    error('Input incorrect length');
end

% Set up positions
x0 = sim.setup.x0(1);
x1 = x0;
midpoint = (sim.setup.x0 + sim.setup.x1)/2;

% Iterate through each section to build up material
poly = zeros(sim.setup.pdeg1, sim.setup.nx);
for i = 1:linput
    
    if linput > 1
        x1 = x1 + sim.setup.sec(i);
        sec = (midpoint <= x1).*(midpoint >= x0);
        x0 = x1;
    else
        sec = ones(length(sim.setup.x0),1);
    end
    % If values then
    if isa(input(i),'double')
        
        if input(i) == -Inf
            % Load the data file again and build the correct function for
            % the interpolation
            material = sim.material.material;
            material = strsplit(material,'&');
            secmat = strtrim(material{i});
            
            if exist(secmat, 'file') == 2
                loadedMat = load(secmat);
            elseif exist(strcat(secmat,'.mat'), 'file') == 2
                loadedMat = load(strcat(secmat,'.mat'));
            else
                error('Material Does Not Exist');
            end
            
            VEgreq = strsplit(sim.material.VEg,'&');
            VEgreq = VEgreq{i};
            VEgreq = VEgreq(5:end);
            func = strcat('@(x)interp1([', num2str(loadedMat.VEg),'], [', num2str(loadedMat.(field)),'] ,', VEgreq ,')');
            
            poly(:, sec==1) = func2sbp(eval(func), sim, sec)/scale;
            
        else
            poly(end, sec==1) = input(i)/scale;
        end
        
    elseif isa(input(i),'cell')
        
        if isa(input{i}, 'double')
            poly(end, sec==1) = input{i}/scale;
        else
            func = input{i};
            poly(:, sec==1) = func2sbp(eval(func), sim, sec)/scale;
        end
        
    else
        class(input)
        
        error('Input must be a list of doubles, or an anonymous function string');
    end
    
end

