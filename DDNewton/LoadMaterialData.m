function sim = LoadMaterialData(sim)

material = sim.material.material;

if ischar(material)
    material = strsplit(material,'&');
    for sec=1:length(material)
        secmat = strtrim(material{sec});
        
        if strcmpi(secmat,'user')
            % Leave material as it is defined on input
        elseif exist(secmat, 'file') == 2
            loadedMat = load(secmat);
            sim.material = buildMat(sim.material, loadedMat, sec);
        elseif exist(strcat(secmat,'.mat'), 'file') == 2
            loadedMat = load(strcat(secmat,'.mat'));
            sim.material = buildMat(sim.material, loadedMat, sec);
        else
            error('Material Does Not Exist');
        end
    end
end

end

function mat = buildMat(mat, mat_in, sec)
% We have a loaded material mat_in which we want to make mat

Eg_req = mat.VEg0;

if ~ischar(Eg_req) % Create Sections
    if isfield(mat_in, 'VEg')
        if length(mat_in.VEg) > 1 % Interpolate
            mat = interpData(Eg_req(sec), mat, mat_in, sec);
        else % Copy
            mat = copyData(mat, mat_in, sec);
        end
    end
else % Set params to -inf and attach data for later conversion to sbp
    
    len = length(Eg_req);
    Eg_req = strsplit(Eg_req,'&');
    
    if len~= length(Eg_req)
        Eg_req = strtrim(Eg_req{sec});
    end
    
    if ~isnan(str2double(Eg_req)) % If number as string then treat as number
        Eg_req = str2double(Eg_req);
        if length(mat_in.VEg) > 1 % Interpolate
            mat = interpData(Eg_req, mat, mat_in, sec);
        else % Copy
            mat = copyData(mat, mat_in, sec);
        end
        
    else % We have a function
        f = fields(mat_in);
        for i = 1:length(f)
            variable = f{i};      
            if ~strcmp(variable, 'VEg')
                if ~isequal(repelem(mat_in.(variable)(1), length(mat_in.(variable))), mat_in.(variable))
                    mat.(variable)(sec)= -inf;
                end
            end
        end
    end
    
end
end

function mat = interpData(Eg_req, mat, mat_in, sec)

Eg_in = mat_in.VEg;
f = fields(mat_in);

for i = 1:length(f)
    if ~strcmp(f{i}, 'VEg')
        variable = mat_in.(f{i});
        if ~strcmp(f{i}, 'eps') && ~strcmp(f{i}, 'reference')
            mat.(f{i})(sec)= interp1(Eg_in, variable, Eg_req);
        end
    end
end
end

function mat = copyData(mat, mat_in, sec)
f = fields(mat_in);
for i = 1:length(f)
    if ~strcmp(f{i}, 'VEg') && ~strcmp(f{i}, 'eps')  && ~strcmp(f{i}, 'reference')
        if length(mat_in.(f{i})) == 1
            mat.(f{i})(sec) = mat_in.(f{i});
        end
    end
end
end