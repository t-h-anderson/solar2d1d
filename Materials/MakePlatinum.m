
clear();

platinum = importdata('./Materials/DataFiles/Pt_Werner_2009.csv');

eps = [platinum.data(:,1), (platinum.data(:,2) + 1i*platinum.data(:,3)).^2];
reference = 'W. S. M. Werner, K. Glantschnig, C. Ambrosch-Draxl. Optical constants and inelastic electron-scattering data for 17 elemental metals, J. Phys Chem Ref. Data 38, 1013-1092 (2009)';
lambda0 = min(platinum.data(:,1));
lambda1 = max(platinum.data(:,1));

clear('platinum');

save('./Materials/Metals/Platinum.mat')