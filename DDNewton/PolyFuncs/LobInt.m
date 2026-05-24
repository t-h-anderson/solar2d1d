%%%%%%%%%%%%%%%%%%%%LobtInt%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ integral ] = LobInt(func, vars, sim)
% Integrate a function using legnedre points

wi = sim.local.wi;
      
if size(vars,1) == sim.setup.pdeg1

    integral =  wi .* func(vars); 

elseif size(vars, 1) == sim.setup.np
    
integral = zeros(sim.setup.np, 1);
    for i = 1:sim.setup.nx
        index1 = (i-1)*sim.setup.pdeg1 + 1;
        index2 = index1 + sim.setup.pdeg1 - 1;
        subvars = vars(index1:index2, :);
        integral(index1:index2) = LobInt(func, subvars, sim);
        
        %integral(index1:index2) = sim.setup.poly*integral(index1:index2);
    end
    
    
    
else
    error('LobInt not defined for this lenght of input')
end
    

end