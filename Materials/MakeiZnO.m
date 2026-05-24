clear();

eps = xlsread('./Materials/DataFiles/iZnO.xlsx');
eps = [eps(:,1), (eps(:,2) + 1i*eps(:,3)).^2];

VEg=[3.3];
VChi=[4.4];
icm3Nc = [3e18];
icm3Nv = [1.7e19];
cm2iVismun = [100];
cm2iVismup = [31];
epsdc = [9];
alpha = [1e-10];
istausrhn = [2e-10];
istausrhp = [1e-8];

save('./Materials/Dielectrics/ZnO.mat')