
clear();

n = xlsread('./Materials/DataFiles/AZOold.xlsx');
eps = [n(:,1), (n(:,2)).^2];

save('./Materials/Dielectrics/AZOold.mat')