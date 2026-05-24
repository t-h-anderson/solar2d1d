xin = [0.2,0.5];
yin = 0.9;

Ep = [0.996, 1.281, 1.860, 2.909, 3.050, 3.599, 4.709, 5.216, 6.399;
    1.311, 1.504, 2.769, 3.008, 3.336, 3.662, 4.919, 5.511, 6.588;
    1.445, 1.523, 2.619, 3.083, 3.212, 3.787, 4.962, 5.638, 6.655;
    1.713, 1.855, 3.024, 3.227, 3.582, 3.97, 5.104, 5.759, 6.729;
    1.268, 1.453, 2.449, 3.005, 3.162, 3.696, 4.885, 5.424, 6.896;
    1.355, 1.595, 2.948, 2.986, 3.332, 3.563, 4.95, 5.532, 6.509;
    1.449, 3.059, 3.553, 4.918, 6.02, 0, 0, 0, 0];

A = [12.255, 14.58, 21.014, 13.185, 58.905, 35.077, 12.44, 7.692, 54.873;
    20.119, 17.816, 17.292, 31.278, 4.029, 91.535, 7.514, 7.366, 40.931;
    19.547, 33.222, 21.356, 55.039, 4.694, 80.262, 9.025, 3.579, 51.707;
    13.966, 34.567, 22.058, 41.486, 18.203, 107.55, 14.859, 55.933, 51.6;
    29.034, 30.033, 44.324, 22.159, 3.116, 47.589, 7.36, 42.131, 35.97;
    15.013, 23.516, 18.503, 21.522, 3.42, 148.907, 7.182, 3.674, 48.984;
    22.696, 28.043, 54.974, 23.104, 88.14, 0, 0, 0, 0];

C = [0.175, 1.087, 2.446, 0.802, 0.425, 1.332, 1.129, 0.928, 3.381;
    0.292, 0.636, 2.088, 0.677, 0.86, 0.901, 1.019, 0.727, 3.189;
    0.197, 0.813, 1.695, 0.871, 0.648, 1.087, 1.092, 0.754, 3.16;
    0.252, 1.162, 0.946, 0.487, 0.733, 0.955, 1.286, 0.592, 3.063;
    0.302, 0.583, 2.837, 0.676, 1.127, 0.67, 1.006, 0.702, 3.772;
    0.612, 1.791, 1.439, 0.591, 0.687, 1.126, 1.102, 0.775, 3.519;
    1.305, 1.269, 1.389, 2.081, 2.956, 0, 0, 0, 0];

Eg = [0.947, 1.001, 1.438, 1.485, 3.017, 2.593, 2.013, 2.324, 2.948;
    1.254, 1.43, 1.366, 2.272, 1.746, 3.302, 1.436, 2.97, 2.367;
    1.397, 1.5, 1.801, 2.395, 2.356, 3.274, 1.637, 1.708, 2.939;
    1.626, 1.728, 2.098, 2.678, 2.659, 3.441, 1.946, 4.973, 3.245;
    1.237, 1.429, 1.647, 2.136, 1.247, 3.312, 1.444, 4.46, 2.286;
    1.287, 1.453, 1.901, 2.395, 1.428, 3.284, 1.415, 2.477, 2.553;
    1.349, 2.019, 2.986, 1.455, 4.315, 0, 0, 0, 0];

% Create function for each sample
eps2list={};
for i = 1:7
    eps2list{i} = @(E) 0;
    for j = 1:length(Eg(i,:))
        eps2list{i} =@(E)( eps2list{i}(E) + (E > Eg(i,j)) .* A(i,j) .* C(i,j).*Ep(i,j).*(E-Eg(i,j)).^2./((E.^2-Ep(i,j).^2).^2.*E+C(i,j).^2 *E.^3));
    end
end

% Find range for x input
x = [0,0.4, 0.63, 1];
j0 = @(E) 4 - (E<0.4) - (E<0.63) - (E<1);
j1 = @(E) 5 - (E<0.4) - (E<0.63) - (E<=1);

eps2lower = eps2list(j0(xin));
eps2upper = eps2list(j1(xin));

xi = @(xin) (j0(xin) ~= j1(xin)) .*  (xin - x(j0(xin)))./(x(j1(xin)) - x(j0(xin)));


% Set peaks correctly
E0 = 1 + 0.71 * xin + 0.34 * (0.9-yin);
E1A = 2.94 + 0.39 * xin;
E1B = 3.71 + 0.49 * xin;
E2B = 4.71 + 0.44 * xin;
E3 = 5.24 + 0.64 * xin;

delta_E0lower = abs(E0 - Ep(j0(xin), 1).');
delta_E0upper = abs(E0 - Ep(j1(xin), 1).');

delta_E1Alower = abs(E1A - Ep(j0(xin), 4).');
delta_E1Aupper = abs(E1A - Ep(j1(xin), 4).');

delta_E1Blower = abs(E1B - Ep(j0(xin), 6).');
delta_E1Bupper = abs(E1B - Ep(j1(xin), 6).');

delta_E2Blower = abs(E2B - Ep(j0(xin), 7).');
delta_E2Bupper = abs(E2B - Ep(j1(xin), 7).');

delta_E3lower = abs(E3 - Ep(j0(xin), 8).');
delta_E3upper = abs(E3 - Ep(j1(xin), 8).');

deltalower =@(E) E + (E<E0).*delta_E0lower +... 
    (E>=E0).*(E<E1A).*(delta_E0lower.*(E - E1A)./(E0 - E1A) + delta_E1Alower.*(E-E0)./(E1A-E0))+ ...
    (E>=E1A).*(E<E1B).*(delta_E1Alower.*(E-E1B)./(E1A - E1B) + delta_E1Blower.*(E-E1A)./(E1B-E1A))+...
    (E>=E1B).*(E<E2B).*(delta_E1Blower.*(E-E2B)./(E1B - E2B) + delta_E2Blower.*(E-E1B)./(E2B-E1B))+...
    (E>=E2B).*(E<E3).*(delta_E2Blower.*(E-E3)./(E2B - E3) + delta_E3lower.*(E-E2B)./(E3-E2B))+...
    (E>=E3).*(delta_E3lower);

deltaupper =@(E) E - (E<E0).*delta_E0upper +... 
    (E>=E0).*(E<E1A).*(delta_E0upper.*(E - E1A)./(E0 - E1A) + delta_E1Aupper.*(E-E0)./(E1A-E0))+ ...
    (E>=E1A).*(E<E1B).*(delta_E1Aupper.*(E-E1B)./(E1A - E1B) + delta_E1Bupper.*(E-E1A)./(E1B-E1A))+...
    (E>=E1B).*(E<E2B).*(delta_E1Bupper.*(E-E2B)./(E1B - E2B) + delta_E2Bupper.*(E-E1B)./(E2B-E1B))+...
    (E>=E2B).*(E<E3).*(delta_E2Bupper.*(E-E3)./(E2B - E3) + delta_E3upper.*(E-E2B)./(E3-E2B))+...
    (E>=E3).*(delta_E3upper);


% eps2l = @(E) cellfun(@(f)f(deltalower(E)),eps2lower, 'UniformOutput',false);
 %eps2u = @(E) cellfun(@(f)f(deltaupper(E)),eps2upper, 'UniformOutput',false);
 
 for i = 1:length(xin)
     fnl = eps2lower{i};
     eps2{i} = fnl(deltalower(ee'))
 end

eps2 = @(E) cell2mat(eps2l(E)) .* (1 - xi(xin)) + xi(xin) .* cell2mat(eps2u(E));

faiz=xlsread('CIGSxt0p2.xlsx');

eps2f = @(E) interp1(faiz(:,1), faiz(:,3), E);

disp();