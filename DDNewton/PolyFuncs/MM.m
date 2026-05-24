function [M] = MM(w,p,q,sim)
% p is a list of polynomials
% q is a list of polynomials
% Creates the matrix of \int_-1^1 p_i q_j

ip = size(p,2);
jq = size(q,2);

M = zeros(ip,jq);


wsum = sum(abs(w));


mult = 1;
if wsum == 0
    return
elseif wsum == w(end)
    q = q*w(end);
    mult = 0;
end

if isequal(p,q)
    for i = 1:ip
        for j = i:jq
            
            p0 = p(:,i);
            
            if mult == 1
                q0 = polytimes(w, q(:,j));
            else
                q0 = q(:,j);
            end
            
            M(i,j) = pint3(p0, q0, sim);
            M(j,i) = M(i,j);
        end
    end
    
else
    for i = 1:ip
        for j = 1:jq
            p0 = p(:,i);
            
            if mult == 1
                q0 = polytimes(w, q(:,j));
            else
                q0 = q(:,j);
            end
            
            M(i,j) = pint3(p0, q0, sim);
        end
    end
end

end
