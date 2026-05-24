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
            out(:,i) = padarray(temp,[0,pdeg-size(temp,2)],0,'pre');
        end
    end
end

end

