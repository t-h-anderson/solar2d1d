clear();

reference = 'Optics, Minoura 2015, https://doi.org/10.1063/1.4921300, Electrical, Bhuiyan 2012 CIGS modelling, Markus Thesis';
VEg=[1.13];
VChi=[4.4];
icm3Nc = [7e18];
icm3Nv = [4.5e18];
cm2iVismun = [40];
cm2iVismup = [12.6];
epsdc = [14];
eps = '@(nmlambda, VEg) CIGs(nmlambda, VEg)';
alpha = [1e-10];
istausrhn = [1e-14];
istausrhp = [1e-14];

save('./Materials/Semiconductors/CIGS1.mat')