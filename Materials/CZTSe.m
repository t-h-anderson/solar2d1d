close all;
clc;

eps11=[];
eps22=[];
 for xin=0:0.1:1   
lambda=[300:10:1400];
CZTS=xlsread(fullfile(fileparts(mfilename('fullpath')), 'DataFiles', 'CZTS.xlsx'));
ReCZTS=interp1(CZTS(:,1),CZTS(:,2),lambda);
ImCZTS=interp1(CZTS(:,1),CZTS(:,3),lambda);

Ep = [1.177, 2.401, 3.789, 4.549,5.886].';

A = [11.12, 25.14, 63.40, 173.08, 231.27].';

C = [1.33, 1.664, 1.563, 2.205, 2.009].';

Eg = [0.7, 1.15,2.598, 3.566, 5.455].';

epsinf = 0.77;
eelower=1240./lambda-0.34*xin;
gamma = sqrt(Ep.^2 - C.^2/2);
beta = sqrt(4*Ep.^2 - C.^2);
eps1lower=zeros(1,length(eelower));
eps2lower=zeros(1,length(eelower));
test1lower=zeros(1,length(eelower));
test2lower=zeros(1,length(eelower));
test3lower=zeros(1,length(eelower));
test4lower=zeros(1,length(eelower));
for i = 1:5
   eps2lower =eps2lower+ (eelower > Eg(i)) .*A(i) .* C(i).*Ep(i).*(eelower - Eg(i)).^2./((eelower.^2-Ep(i).^2).^2.*eelower + C(i).^2 .*eelower.^3); 
end

for i = 1:5
   xi4lower = (eelower.^2 - gamma(i).^2).^2 + beta(i).^2 .* C(i).^2/4;
  
   atanlower = (eelower.^2 - Ep(i).^2) .* (Ep(i).^2 + Eg(i).^2) + Eg(i).^2 .* C(i).^2;
  
   
   alnlower = (Eg(i).^2 - Ep(i).^2) .* eelower.^2 + Eg(i).^2 .* C(i).^2 - Ep(i).^2 .*(Ep(i).^2 + 3 * Eg(i).^2);
  
   
   eps1lower = eps1lower +  ...
       (A(i).*C(i) .* alnlower./(2 * pi *  xi4lower .* beta(i) .* Ep(i))) .* log((Ep(i).^2 +Eg(i).^2 + beta(i) .* Eg(i))./((Ep(i).^2 +Eg(i).^2 - beta(i) .* Eg(i))))-...
      (A(i) .* atanlower./(pi *  xi4lower .* Ep(i))) .* (pi - atan((2*Eg(i) + beta(i))./C(i)) + atan((- 2*Eg(i) + beta(i))./C(i))   )+...
     2*(A(i) .* Ep(i) .* Eg(i) .* (eelower.^2 - gamma(i).^2))./(pi *  xi4lower .* beta(i)) .* (pi + 2* atan(2*(gamma(i).^2 - Eg(i).^2)./(beta(i).*C(i)))  )-...
    (A(i) .* Ep(i) .* C(i) .* (eelower.^2 + Eg(i).^2))./(pi *  xi4lower .* eelower) .* log(abs(eelower-Eg(i))./(eelower + Eg(i)))+...
     2*(A(i) .* Ep(i) .* C(i) .* Eg(i))./(pi *  xi4lower) .* log(abs(eelower-Eg(i)).*(eelower + Eg(i))./sqrt((Ep(i).^2-  Eg(i).^2).^2 + Eg(i).^2.*C(i).^2)); 
    

end
eps1lower =eps1lower+epsinf ;
eps1 = eps1lower *(1-xin) + ReCZTS *(xin);
eps2 = eps2lower *(1-xin) + ImCZTS *(xin);
%plot(eps1)
%hold on
%plot(eps2)
%hold on
eps11=[eps11;eps1];
eps22=[eps22;eps2];
 end
csvwrite('eps1.csv',eps11.');

csvwrite('eps2.csv',eps22.');
