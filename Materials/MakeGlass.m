
clear();

glass = importdata('./Materials/DataFiles/N-BK7.csv');

n = [glass.data(:,1), glass.data(:,2)];
k = [glass.data(:,3), glass.data(:,4)];
    
k = [glass.data(:,1), interp1(k(~isnan(k(:,1)),1), k(~isnan(k(:,1)),2), n(:,1))];

eps = [glass.data(:,1), (n(:,2)+k(:,2)*1i).^2];

reference = 'refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT';
lambda0 = min(eps(:,1));
lambda1 = max(eps(:,1));

clear('glass');

save('./Materials/Dielectrics/Glass.mat')