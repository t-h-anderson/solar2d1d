clear();
reference = 'Optical: J. Phys.:Conf. Ser. 286, 012038 (2011), Electrical: Frisk et al, Appl. Phys. 47, 485104 (2014)';
eps = xlsread('./Materials/DataFiles/CdS.xlsx');

VEg = [2.4];
VChi=[4.2];
icm3Nc = [1.3e18];
icm3Nv = [9.1e19];
cm2iVismun = [72];
cm2iVismup = [20];
epsdc = [5.4];
eps = [eps(:,1), (eps(:,2) + 1i*eps(:,3)).^2];
alpha = [0e-10];
istausrhn = [2e-9];
istausrhp = [1e-6];

save('./Materials/Semiconductors/CdS.mat')