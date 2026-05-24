function [P] = lobpoly(n)
% Creates matrix with columns given by n lobato polynomials

[xi,~]=lgnodes(n);

xi = fliplr(xi);

P=zeros(n+1,n+1);

for i=1:n+1
    P(:,i) = (poly([xi(1:i-1),xi(i+1:end)]).');
    P(:,i) = P(:,i)./polyval(P(:,i),xi(i));
end
end

