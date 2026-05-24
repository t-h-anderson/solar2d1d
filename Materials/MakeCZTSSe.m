clear();

reference = 'Optics, Minoura 2015, https://doi.org/10.1063/1.4921300, Electrical, Bhuiyan 2012 CIGS modelling, Markus Thesis';

VEg=[1.0];
VChi=[0];
icm3Nc = [6.8e17];
icm3Nv = [1.5e19];
cm2iVismun = [100];
cm2iVismup = [13];
epsdc = [13.6];
eps = '@(nmlambda, VEg) CZTSSe(nmlambda, VEg)';
alpha = [0e-10];
istausrhn = [1e-9];
istausrhp = [1e-6];

save('./Materials/Semiconductors/CZTSSe.mat')