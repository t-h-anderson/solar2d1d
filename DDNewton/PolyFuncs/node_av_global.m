function [ node_vals ] = node_av_global(p)
% Take a standard basis piecewise polynomial, and produce average values
% for the endpoints

nx = size(p,2);

node_vals = zeros(nx-1,1);

for i = 1:nx-1
    node_vals(i) = 0.5*(polyval(p(:,i),1) + polyval(p(:,i+1),-1));
end
    

end

