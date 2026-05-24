function [ answer] = tiprods(X, Y, sim)
% Calculate the inner product of a rank 3 tensor whic h is block diagonal
% and is symmetric with each index

np = sim.setup.np;
pdeg1 = sim.setup.pdeg1;
nx = sim.setup.nx;

% Take Y and make it repeat pdeg1 times every pdeg1 elements, e.g.
% [1,2,3,4]->[1,2,1,2,3,4,3,4] if pdeg1 = 2
Ytemp=coeffs2ibp(Y, sim);
Y2 = reshape(repmat(Ytemp,[pdeg1,1]),[np*pdeg1,1]);

% Form a diagonal sparse matrix where each column is an pdeg1 block of Y,
% the next pdeg columns are the same but pdeg1 lower. e.g. if Y = [1,2,3,4]
% then V2 = [1 0 0 0; 2 0 0 0; 0 1 0 0; 0 2 0 0; 0 0 3 0; 0 0 4 0; 0 0 0 3; 0 0 0 4] 
index = 1:np*pdeg1;
col = reshape(kron(1:nx*pdeg1,ones(pdeg1,1)),[np*pdeg1,1]);
V2 = sparse(index, col, Y2);

answer = X*V2;


end

