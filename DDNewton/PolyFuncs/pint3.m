function I = pint3(p, q, sim)

    np = size(p,1);
    nq = size(q,1);
    
    n = sim.setup.pdeg;
    
    
    toggle = 0;
    switch (np-1) 
        case n
            toggle = toggle + 1;
        case 2*n
            toggle = toggle + 2;
    end
    
    switch (nq-1) 
        case n
            toggle = toggle + 10;
        case 2*n
            toggle = toggle + 20;
    end
    
    %disp('Pint3 set to slow');
    %toggle = 0;
    
    switch toggle
        case 11
            mat = sim.int.nnmat;
        case 12
            mat = sim.int.n2nmat';
        case 21
            mat = sim.int.n2nmat;
        case 22
            mat = sim.int.n2n2mat;
        otherwise
            mat = zeros(np, nq);
            vect = zeros(1,2*max(np,nq) - 1);  
            vect(1:2:2*max(np,nq)-1) = 1./(2*max(np,nq)-1:-2:0);
            for i=1:nq
                mat(:,i) = vect(end - (nq) + i - (np) + 1 : end - (nq) + i);
            end
    end
    
    
    I = 2*p'*mat*q;    
    
end
