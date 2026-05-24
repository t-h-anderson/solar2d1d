function [ answer] = tiprod(X, Xind, Y)
% Calculate the inner product of a rank 3 tensor with a vector, with Xind the
% index to sum over.

s = size(X);

if length(s)~=3
   error('Tensor must be rank 3');
end

switch Xind
    case 1
        answer = zeros(s(2),s(3));
        
        if length(Y)~= s(1)
            error('Tensor and vector are not compatible sizes in this dimension')
        end
        
        for i = 1:s(1)
            answer = answer + squeeze(X(i,:,:))*Y(i);
        end
    case 2
        answer = zeros(s(1),s(3));
        
        if length(Y)~= s(2)
            error('Tensor and vector are not compatible sizes in this dimension')
        end
        
        for i = 1:s(2)
            answer = answer + squeeze(X(:,i,:))*Y(i);
        end
    case 3
        answer = zeros(s(1),s(2));
        
        if length(Y)~= s(3)
            error('Tensor and vector are not compatible sizes in this dimension')
        end
        
        
        for i = 1:s(3)
            answer = answer + squeeze(X(:,:,i))*Y(i);
        end
    otherwise
        error('Requested dimension must be 1, 2, or 3.')
end
end

