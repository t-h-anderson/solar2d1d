function [ roots ] = polyrootcount(p, val)
% Count the number of roots in a piecewise polynomial p - val
% between -1-delta and 1+delta

p(end,:) = p(end,:) - val;
polys{1} = p;
polys{2} = pdiff(p, 1, 0);

nextp = inf;
i = 1;
while sum(sum(nextp)) ~= 0
    
    nextp = 0*polys{1};
    
    % Calculate the Strum sequence for each part of the piecewise polynomial
    p1 = polys{i};
    p2 = polys{i+1};    
    for j = 1:size(p,2)
        if sum(p2(:,j)) ~=0
            [~, nextp(:,j)] = deconv(p1(:,j), p2(find(p2(:,j)~=0,1):end,j));
        end
    end
    polys{i+2} = - nextp;
    
    i = i + 1;
end

for j = 1:i
    poly = polys{j};
    for k = 1:size(poly, 2)
        a(j, k) = polyval(poly(:,k), -1 - 1e-9);
        b(j, k) = polyval(poly(:,k), 1.0 + 1e-9);
    end
end

for k = 1:size(poly, 2)
    aval(k) = sum(abs(conv(sign(a(:,k))',[0.5;-0.5],'valid')));
    bval(k) = sum(abs(conv(sign(b(:,k))',[0.5;-0.5],'valid')));
end

roots = abs(aval - bval);

for k = 1:size(poly, 2)
    if roots(k) == 0 && a(1,k) < 0
        roots(k) = inf;
    end
end


end

