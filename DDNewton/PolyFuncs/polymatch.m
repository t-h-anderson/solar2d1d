function [ p ] = polymatch(p, n)
% Make an polynomial array p an array of poly of degree n

if size(p,1)<n
    p = padarray(p, [n - size(p,1), 0], 0, 'pre');
else
    p = p(end-n+1:end, :);
end

end

