
clear();

eps = xlsread('./Materials/DataFiles/nAZO.xlsx');
eps = [eps(:,1), (eps(:,2)).^2];

save('./Materials/Dielectrics/AZOnoabs.mat')