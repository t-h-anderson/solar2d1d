function [ p ] = polymatch(p, n)
% Make an polynomial array p an array of poly of degree n

if size(p,1)<n
    % padarray(p, [n - size(p,1), 0], 0, 'pre'): prepend zero rows so the
    % polynomial has degree n-1. Inlined to avoid the Image Processing
    % Toolbox dependency.
    p = [zeros(n - size(p,1), size(p,2)); p];
else
    p = p(end-n+1:end, :);
end

end

