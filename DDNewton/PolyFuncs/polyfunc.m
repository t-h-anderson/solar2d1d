function [ pout ] = polyfunc(fn, p)
% polyfunc(fn, pin, h) applies a function 'fn' to the piecewise
% polynomials 'p'(xi), where xi in (-1,1)
% where 'h' is the vector of mesh sizes.


[pdeg1, nx] = size(p);
res = (pdeg1);

xbar = lgnodes(res-1);
%xbar = linspace(-1,1,res+10);


pout = zeros(pdeg1, nx);
for i = 1:nx
    y = fn(polyval(p(:,i), xbar));
    pout(:, i) = polyfit_nowarn(xbar, y, pdeg1-1);
end

end

