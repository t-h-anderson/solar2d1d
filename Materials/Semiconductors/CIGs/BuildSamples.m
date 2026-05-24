function eps = BuildSamples(lambda, VEg)

if size(lambda,2)>size(lambda,1)
    lambda = lambda.';
end

if size(VEg,1)>size(VEg,2)
    VEg = VEg.';
end

yin = 0.9;

ee = 1240./lambda;
xin = interp1([0.947, 1.254, 1.397, 1.626], [0, 0.4, 0.63, 1], VEg);

Ep = [0.996, 1.281, 1.860, 2.909, 3.050, 3.599, 4.709, 5.216, 6.399;
    1.311, 1.504, 2.769, 3.008, 3.336, 3.662, 4.919, 5.511, 6.588;
    1.445, 1.523, 2.619, 3.083, 3.212, 3.787, 4.962, 5.638, 6.655;
    1.713, 1.855, 3.024, 3.227, 3.582, 3.97, 5.104, 5.759, 6.729;
    1.268, 1.453, 2.449, 3.005, 3.162, 3.696, 4.885, 5.424, 6.896;
    1.355, 1.595, 2.948, 2.986, 3.332, 3.563, 4.95, 5.532, 6.509;
    1.449, 3.059, 3.553, 4.918, 6.02, 0, 0, 0, 0].';

A = [12.255, 14.58, 21.014, 13.185, 58.905, 35.077, 12.44, 7.692, 54.873;
    20.119, 17.816, 17.292, 31.278, 4.029, 91.535, 7.514, 7.366, 40.931;
    19.547, 33.222, 21.356, 55.039, 4.694, 80.262, 9.025, 3.579, 51.707;
    13.966, 34.567, 22.058, 41.486, 18.203, 107.55, 14.859, 55.933, 51.6;
    29.034, 30.033, 44.324, 22.159, 3.116, 47.589, 7.36, 42.131, 35.97;
    15.013, 23.516, 18.503, 21.522, 3.42, 148.907, 7.182, 3.674, 48.984;
    22.696, 28.043, 54.974, 23.104, 88.14, 0, 0, 0, 0].';

C = [0.175, 1.087, 2.446, 0.802, 0.425, 1.332, 1.129, 0.928, 3.381;
    0.292, 0.636, 2.088, 0.677, 0.86, 0.901, 1.019, 0.727, 3.189;
    0.197, 0.813, 1.695, 0.871, 0.648, 1.087, 1.092, 0.754, 3.16;
    0.252, 1.162, 0.946, 0.487, 0.733, 0.955, 1.286, 0.592, 3.063;
    0.302, 0.583, 2.837, 0.676, 1.127, 0.67, 1.006, 0.702, 3.772;
    0.612, 1.791, 1.439, 0.591, 0.687, 1.126, 1.102, 0.775, 3.519;
    1.305, 1.269, 1.389, 2.081, 2.956, 0, 0, 0, 0].';

Eg = [0.947, 1.001, 1.438, 1.485, 3.017, 2.593, 2.013, 2.324, 2.948;
    1.254, 1.43, 1.366, 2.272, 1.746, 3.302, 1.436, 2.97, 2.367;
    1.397, 1.5, 1.801, 2.395, 2.356, 3.274, 1.637, 1.708, 2.939;
    1.626, 1.728, 2.098, 2.678, 2.659, 3.441, 1.946, 4.973, 3.245;
    1.237, 1.429, 1.647, 2.136, 1.247, 3.312, 1.444, 4.46, 2.286;
    1.287, 1.453, 1.901, 2.395, 1.428, 3.284, 1.415, 2.477, 2.553;
    1.349, 2.019, 2.986, 1.455, 4.315, 0, 0, 0, 0].';

epsinf = [1.342, 1.374, 1.302, 1.202, 1.236, 1.203, 1.288];

gamma = sqrt(Ep.^2 - C.^2/2);
beta = sqrt(4*Ep.^2 - C.^2);

% Set peaks correctly
E0 = 1 + 0.71 * xin + 0.34 * (0.9-yin);
E1A = 2.94 + 0.39 * xin;
E1B = 3.71 + 0.49 * xin;
E2B = 4.71 + 0.44 * xin;
E3 = 5.24 + 0.64 * xin;

% Find two samples required for x input
x = [0,0.4, 0.63, 1];
j0 = 4 - (xin<0.4) - (xin<0.63) - (xin<1);
j1 = 5 - (xin<0.4) - (xin<0.63) - (xin<=1);
xi = (j0 ~= j1) .*  (xin - x(j0))./(x(j1) - x(j0));
xi(isnan(xi))= 1;

% Create deltas for both samples
delta_E0lower = (E0 - Ep(1, j0));
delta_E0upper = (E0 - Ep(1, j1));

delta_E1Alower = (E1A - Ep(4, j0));
delta_E1Aupper = (E1A - Ep(4, j1));

delta_E1Blower = (E1B - Ep(6, j0));
delta_E1Bupper = (E1B - Ep(6, j1));

delta_E2Blower = (E2B - Ep(7, j0));
delta_E2Bupper = (E2B - Ep(7, j1));

delta_E3lower = (E3 - Ep(8, j0));
delta_E3upper = (E3 - Ep(8, j1));

eelower = ee - ((ee<E0).*delta_E0lower +... 
    (ee>=E0).*(ee<E1A).*(delta_E0lower.*(ee - E1A)./(E0 - E1A) + delta_E1Alower.*(ee-E0)./(E1A-E0))+ ...
    (ee>=E1A).*(ee<E1B).*(delta_E1Alower.*(ee-E1B)./(E1A - E1B) + delta_E1Blower.*(ee-E1A)./(E1B-E1A))+...
    (ee>=E1B).*(ee<E2B).*(delta_E1Blower.*(ee-E2B)./(E1B - E2B) + delta_E2Blower.*(ee-E1B)./(E2B-E1B))+...
    (ee>=E2B).*(ee<E3).*(delta_E2Blower.*(ee-E3)./(E2B - E3) + delta_E3lower.*(ee-E2B)./(E3-E2B))+...
    (ee>=E3).*(delta_E3lower));

eeupper = ee - ((ee<E0).*delta_E0upper +... 
    (ee>=E0).*(ee<E1A).*(delta_E0upper.*(ee - E1A)./(E0 - E1A) + delta_E1Aupper.*(ee-E0)./(E1A-E0))+ ...
    (ee>=E1A).*(ee<E1B).*(delta_E1Aupper.*(ee-E1B)./(E1A - E1B) + delta_E1Bupper.*(ee-E1A)./(E1B-E1A))+...
    (ee>=E1B).*(ee<E2B).*(delta_E1Bupper.*(ee-E2B)./(E1B - E2B) + delta_E2Bupper.*(ee-E1B)./(E2B-E1B))+...
    (ee>=E2B).*(ee<E3).*(delta_E2Bupper.*(ee-E3)./(E2B - E3) + delta_E3upper.*(ee-E2B)./(E3-E2B))+...
    (ee>=E3).*(delta_E3upper));


% Create function for each sample
eps1lower = ones(length(ee), length(xin));
eps1upper = ones(length(ee), length(xin));

test1lower = zeros(length(ee), length(xin));
test1upper = zeros(length(ee), length(xin));

test2lower = zeros(length(ee), length(xin));
test2upper = zeros(length(ee), length(xin));

test3lower = zeros(length(ee), length(xin));
test3upper = zeros(length(ee), length(xin));

test4lower = ones(length(ee), length(xin));
test4upper = ones(length(ee), length(xin));

for i = 1:9
   xi4lower = (eelower.^2 - gamma(i,j0).^2).^2 + beta(i,j0).^2 .* C(i,j0).^2/4;
   xi4upper = (eeupper.^2 - gamma(i,j1).^2).^2 + beta(i,j1).^2 .* C(i,j1).^2/4; 
  
   atanlower = (eelower.^2 - Ep(i,j0).^2) .* (Ep(i,j0).^2 + Eg(i,j0).^2) + Eg(i,j0).^2 .* C(i,j0).^2;
   atanupper = (eeupper.^2 - Ep(i,j1).^2) .* (Ep(i,j1).^2 + Eg(i,j1).^2) + Eg(i,j1).^2 .* C(i,j1).^2;
   
   alnlower = (Eg(i,j0).^2 - Ep(i,j0).^2) .* eelower.^2 + Eg(i,j0).^2 .* C(i,j0).^2 - Ep(i,j0).^2 .*(Ep(i,j0).^2 + 3 * Eg(i,j0).^2);
   alnupper = (Eg(i,j1).^2 - Ep(i,j1).^2) .* eeupper.^2 + Eg(i,j1).^2 .* C(i,j1).^2 - Ep(i,j1).^2 .*(Ep(i,j1).^2 + 3 * Eg(i,j1).^2);
   
   eps1lower = eps1lower +  ...
       (A(i,j0).*C(i,j0) .* alnlower./(2 * pi *  xi4lower .* beta(i,j0) .* Ep(i,j0))) .* log((Ep(i,j0).^2 +Eg(i,j0).^2 + beta(i,j0) .* Eg(i,j0))./((Ep(i,j0).^2 +Eg(i,j0).^2 - beta(i,j0) .* Eg(i,j0))));
    test1lower = test1lower - (A(i,j0) .* atanlower./(pi *  xi4lower .* Ep(i,j0))) .* (pi - atan((2*Eg(i,j0) + beta(i,j0))./C(i,j0)) + atan((- 2*Eg(i,j0) + beta(i,j0))./C(i,j0))   );
    test2lower = test2lower + 2*(A(i,j0) .* Ep(i,j0) .* Eg(i,j0) .* (eelower.^2 - gamma(i,j0).^2))./(pi *  xi4lower .* beta(i,j0)) .* (pi + 2* atan(2*(gamma(i,j0).^2 - Eg(i,j0).^2)./(beta(i,j0).*C(i,j0)))  );
    test3lower = test3lower - (A(i,j0) .* Ep(i,j0) .* C(i,j0) .* (eelower.^2 + Eg(i,j0).^2))./(pi *  xi4lower .* eelower) .* log(abs(eelower-Eg(i,j0))./(eelower + Eg(i,j0)));
    test4lower = test4lower + 2*(A(i,j0) .* Ep(i,j0) .* C(i,j0) .* Eg(i,j0))./(pi *  xi4lower) .* log(abs(eelower-Eg(i,j0)).*(eelower + Eg(i,j0))./sqrt((Ep(i,j0).^2-  Eg(i,j0).^2).^2 + Eg(i,j0).^2.*C(i,j0).^2)); 

   eps1upper= eps1upper + ...
       (A(i,j1).*C(i,j1) .* alnupper./(2 * pi *  xi4upper .* beta(i,j1) .* Ep(i,j1))) .* log((Ep(i,j1).^2 +Eg(i,j1).^2 + beta(i,j1) .* Eg(i,j1))./((Ep(i,j1).^2 +Eg(i,j1).^2 - beta(i,j1) .* Eg(i,j1))));
    test1upper = test1upper - (A(i,j1) .* atanupper./(pi *  xi4upper .* Ep(i,j1))) .* (pi - atan((2*Eg(i,j1) + beta(i,j1))./C(i,j1)) + atan((- 2*Eg(i,j1) + beta(i,j1))./C(i,j1))   );
    test2upper = test2upper + 2*(A(i,j1) .* Ep(i,j1) .* Eg(i,j1) .* (eeupper.^2 - gamma(i,j1).^2))./(pi *  xi4upper .* beta(i,j1)) .* (pi + 2* atan(2*(gamma(i,j1).^2 - Eg(i,j1).^2)./(beta(i,j1).*C(i,j1)))  );
    test3upper = test3upper - (A(i,j1) .* Ep(i,j1) .* C(i,j1) .* (eeupper.^2 + Eg(i,j1).^2))./(pi *  xi4upper .* eeupper) .* log(abs(eeupper-Eg(i,j1))./(eeupper + Eg(i,j1)));
    test4upper = test4upper + 2*(A(i,j1) .* Ep(i,j1) .* C(i,j1) .* Eg(i,j1))./(pi *  xi4upper) .* log(abs(eeupper-Eg(i,j1)).*(eeupper + Eg(i,j1))./sqrt((Ep(i,j1).^2- Eg(i,j1).^2).^2 + Eg(i,j1).^2.*C(i,j1).^2)); 


end




eps2lower = zeros(length(ee), length(xin));
eps2upper = zeros(length(ee), length(xin));
for i = 1:9
   eps2lower = eps2lower + (eelower > Eg(i, j0)) .* A(i, j0) .* C(i, j0).*Ep(i, j0).*(eelower - Eg(i, j0)).^2./((eelower.^2-Ep(i, j0).^2).^2.*eelower + C(i, j0).^2 .*eelower.^3);
   eps2upper = eps2upper + (eeupper > Eg(i, j1)) .* A(i, j1) .* C(i, j1).*Ep(i, j1).*(eeupper - Eg(i, j1)).^2./((eeupper.^2-Ep(i, j1).^2).^2.*eeupper + C(i, j1).^2 .*eeupper.^3);
end

eps1 = eps1lower * diag(1-xi) + eps1upper * diag(xi);

eps2 = eps2lower * diag(1-xi) + eps2upper * diag(xi);

%faiz=xlsread('CIGSxt0p2.xlsx'); 
%eps1f = @(E) interp1(faiz(:,1), faiz(:,2), E);

% surf(ee, xin, eps1.')
% xlabel('Photon energy (eV)')
% ylabel('Ga Composition x')
% zlabel('Imag. Rel. Perm.')
% title('Reprodction of Minoura 2015 Figure 5a')
% xlim([0.75,6.5]);

% figure(1)
% plot(ee, eps1upper, ee,test1upper, ee, test2upper, ee, test3upper, ee, test4upper, ee, eps1upper+test1upper+test2upper+test3upper+test4upper, ee, eps1f(ee))
% legend('1st', '2nd', '3rd', '4th', '5th', 'sum')
% 
% figure(2)
% plot(ee, eps1upper+test1upper+test2upper+test3upper+test4upper, ee, eps1f(ee));
% 
% figure(3)
% plot(ee, eps1lower, ee,test1lower, ee, test2lower, ee, test3lower, ee, test4lower, ee, eps1lower+test1lower+test2lower+test3lower+test4lower, ee, eps1f(ee))
% legend('1st', '2nd', '3rd', '4th', '5th', 'sum')
% 
% figure(4)
% plot(ee, eps1lower+test1lower+test2lower+test3lower+test4lower, ee, eps1f(ee))


%plot(ee, eps1, ee, eps1upper, ee, eps1lower, ee, eps1f(ee))
%legend('Tom', 'x=0, E shifted', 'x=0.4, E shifted', 'Faiz', 'Location', 'SouthEast')
%xlim([0.7,6.5]);
eps = eps1 + 1i*eps2;