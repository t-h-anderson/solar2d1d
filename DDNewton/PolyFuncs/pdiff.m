function out = pdiff(p,n, varargin)
% Differentiate a piecewise polynomial 
out = p;

[pdeg,np] = size(p);

padtoggle = 1;
if length(varargin) == 1
    padtoggle = varargin{1};
end

for i = 1:np
    for j=1:n
        temp = polyder(out(:,i));
        if padtoggle
            % padarray(temp, [0, pdeg-size(temp,2)], 0, 'pre'): polyder
            % returns a row shorter than the input; prepend zeros so the
            % derivative slots back into the same fixed-degree column.
            % Inlined to avoid the Image Processing Toolbox.
            out(:,i) = [zeros(1, pdeg-size(temp,2)), temp];
        end
    end
end

end

