
clear();

 % Electrical Data aSiGeC 
VEg = [1.3, 1.803, 1.95];
VChi = [4.17, 4.03, 3.95];
icm3Nc = [2.5e20, 2.5e20, 2.5e20];
icm3Nv = [2.5e20, 2.5e20, 2.5e20];
cm2iVismun = [15, 15, 15];
cm2iVismup = [3, 3, 3];
epsdc = [11.68, 11.68, 11.68];

reference = 'Forouhi - from Faiz';

data = xlsread('./Materials/DataFiles/aSi_Shokeen.xlsx');

nmlambdas=data(:,1);
eps = [nmlambdas, data(:,2) + 1.i*data(:,3)];

lambda0 = min(nmlambdas);
lambda1 = max(nmlambdas);

clear('data');
clear('nmlamdbas');

save('./Materials/Semiconductors/aSi_Shokeen.mat')