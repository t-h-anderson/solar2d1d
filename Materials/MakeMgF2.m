clear();

eps = xlsread('./Materials/DataFiles/MgF2.xlsx');
eps = [eps(:,1), (eps(:,2) + 1i*1e-9).^2];

save('./Materials/Dielectrics/MgF2.mat')