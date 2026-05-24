
clear();

eps = xlsread('./Materials/DataFiles/RefIndexAZO.xlsx');
eps = [eps(:,1), (eps(:,2) + 1i*eps(:,3)).^2];

save('./Materials/Dielectrics/AZO.mat')