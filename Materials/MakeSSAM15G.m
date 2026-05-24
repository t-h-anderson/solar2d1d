clear();
clc;
data = xlsread('./DataFiles/SSAM15G.xlsx');

data = [data(:,1), data(:,2)];

save('./AM15G.mat')