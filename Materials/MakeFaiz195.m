
clear();

 % Electrical Data aSiGeC 
VEg = [1.95];
VChi = [3.95];
icm3Nc = [2.5e20];
icm3Nv = [2.5e20];
cm2iVismun = [15];
cm2iVismup = [3];
epsdc = [11.68];

reference = 'aSi 1.95 - from Faiz';

data = xlsread('./Materials/DataFiles/aSi195_Faiz.xlsx');

nmlambdas=data(:,1);
eps = [nmlambdas, data(:,2) + 1.i*data(:,3)];

lambda0 = min(nmlambdas);
lambda1 = max(nmlambdas);

clear('data');
clear('nmlamdbas');

save('./Materials/Semiconductors/aSi195_Faiz.mat')