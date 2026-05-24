function [ r ] = polytimes(p, q)

if sum(abs(p(1:end-1))) == 0
    r = p(end) * q;
elseif sum(abs(q(1:end-1))) == 0
    r = q(end) * p;
else
    r = conv(p,q);
end

end

