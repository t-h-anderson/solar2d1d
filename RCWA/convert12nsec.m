
%%%%%%%%%%%%%%%%%%%%%%%%RCWA Initialization
function material = convert12nsec(material, nsec)
% Convert a length 1 input to an nsec length
f = fields(material);
% Loop through the fields
for i = 1:length(f)
    % Get the material property
    prop = material.(f{i});
    if ischar(prop) % if we have a character string split at '&'
        propcells = strsplit(prop,'&');
        if length(propcells) == 1 % If length is 1, then duplicate to nsec
            prop = strtrim(propcells{1});
            p2 = prop;
            for j = 1:nsec-1
                prop = strcat(prop, '& ', p2);
            end
            material.(f{i}) = prop;
        elseif length(propcells) ~= nsec
            error('Input "%s" incorrect length', f{i});
        else
            prop = strtrim(propcells{1});
            if ~strcmp(prop(1),'@') && ~isempty(str2num(prop(1)))
                 prop = strcat('@(x)',prop);
            end
            
            for j = 2:nsec
                p2 = strtrim(propcells{j});
                if ~strcmp(p2(1),'@') && ~isempty(str2num(p2(1)))
                    p2 = strcat('@(x)',p2);
                end
                prop = strcat(prop, '& ', p2);
            end
            material.(f{i}) = prop;
        end
    else 
        if length(prop) == 1 && ~strcmp(f{i}, 'nsec') && ~strcmp(f{i}, 'elecsec')% If we have one number, duplicate to nsec (ignore nsec)
            prop = prop*ones(1,nsec);
            material.(f{i}) = prop;
        elseif length(prop) ~= nsec && ~strcmp(f{i}, 'nsec') && ~strcmp(f{i}, 'elecsec')
            error('Input "%s" incorrect length', f{i});
        end
    end
end

end
