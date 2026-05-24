function [M] = Tensor(p,q,r,d,sim)
% p is a list of polynomials
% q is a list of polynomials
% q is a list of polynomials
% Creates the matrix of \int_-1^1 p^(d(1))_i q^(d(2))_j r^(d(3))_k

I = size(p,2);
J = size(q,2);
K = size(r,2);

if length(d) ==3
   p = pdiff(p, d(1));
   q = pdiff(q, d(2));
   r = pdiff(r, d(3));
end
    
M = zeros(I,J,K);
    
for i = 1:I
    p0 = p(:,i);
    for j = 1:J
        q1 = polytimes(p0, q(:,j));            
        for k = 1:K
            M(i,j,k) = pint3(q1, r(:,k), sim);
        end
    end
end

end
