
clear();

eps = load('./Materials/FTOnk1200.csv', '-ascii');
eps = [eps(:,1), (eps(:,2) + 1i*eps(:,3)).^2];

save('ITO.mat')