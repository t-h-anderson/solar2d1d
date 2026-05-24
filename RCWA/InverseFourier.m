
function [ output ] = InverseFourier(x, fcomp, shift, method, varargin)
% Calculate the inverse fourier of a dataset at points x
% x positions lie in range -Lx/2:Lx/2
% Nt terms of fourier sum taken

fcomp =  squeeze(fcomp);

output = zeros(length(x), size(fcomp,1));

for i = 1:size(fcomp,1)
    
    fcomp_temp = fcomp(i,:);
    
    terms = length(fcomp_temp);
    NT = (terms-1)/2;
    
    if (length(varargin)==1)
        temp=floor(varargin{1});
        NT0 = min(NT,temp);
    else
        NT0 = NT;
    end
    
    if(method ==0)
        
        % Change sign of alternate terms
        if mod(NT,2)
            fcomp_temp(2:2:end) = -(fcomp_temp(2:2:end));
        else
            fcomp_temp(1:2:end) = -(fcomp_temp(1:2:end));    
        end
        
        FFTorder = (2*NT+1)*[fcomp_temp(NT+1:NT+1+NT0),fcomp_temp(NT-NT0+1:NT)];
        
        invff = ifft(FFTorder); % Correct Scaling for zeroth order
     
        xin = linspace(x(1),x(end),length(invff)+1);
        output(1:end-1,i) = interpft(invff, length(x)-1);
        output(end,i) = output(1,i);
        
     %   if shift
     %       output(:,i) = circshift(output(:,i),length(x)*shift);
     %   end
   %     output(:,i) = interp1(xin, [invff,invff(1)], x,'nearest');
        
    elseif (method == 1)
        fsum = 0;
        error('not tested');
        for j = -NT0:NT0
            k =  2*pi * (j) /Lx; % add incident angle
            fsum = fsum + fcomp_temp(j+NT+1).*exp(1i*k*x);
        end
        output(:,i) = fsum;
    else
        error('Inverse method not defined');
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

