

function sim  = RCWA(ll, tt, sim)
%  RCWA(inmk0, radtheta, epsf, iepsf, Nt, sim)
% Z is the Z matrix from EMSW
% Tn the tranfter matrix for each mode through each slice
% R the reflected wave
% E is the resulting electric field
% H is the resulting magnetic field
Nt = sim.setup.Nt;
epsf = sim.zfl.eps{ll};
iepsf = sim.zfl.ieps{ll};
radtheta = sim.setup.radtheta(tt);

n0 = sqrt(epsf(1, 2* sim.setup.Nt +1));
inmk0 = n0*sim.setup.inmk0(ll);
nmik0 = 1/inmk0;
nmdz = sim.setup.nmdz;
Nz = length(nmdz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompose light wave vectors (for incident, reflected and transmitted)
% ESW - eqns. 2.93
inmk0x = inmk0*sin(radtheta)+(-Nt:Nt)*2*pi/sim.input.nmLx;
inmk0z = sqrt(inmk0^2-inmk0x.^2);
 
% Create Y matrices depending on polarization (ESW - eqns. 2.128 and 2.129 split into s and p)
% Set eta0 = 1
if (sim.setup.pol == 1) % nonzero right columns of Y matrices
    Ye_inc = - diag(inmk0z/inmk0);
    Ye_ref = - Ye_inc;
    
    Yh_inc = - eye(2*Nt+1);
    Yh_ref = - eye(2*Nt+1);
elseif (sim.setup.pol == 2) % nonzero left columns of Y matrices
    Ye_inc = eye(2*Nt+1); 
    Ye_ref = eye(2*Nt+1);
    
    Yh_inc = - diag(inmk0z/inmk0);
    Yh_ref = - Yh_inc;
else
    error('Set polarization to 1 for p-polarization, and 2 for s-polarization');
end

% Create column Y matrices for incident and transmitted field
Y_inc = [Ye_inc; Yh_inc]; 
Y_ref = [Ye_ref; Yh_ref];

% Create Z matrices
Z  = cell(Nz+1,1);
Z{Nz+1} = Y_inc; %ESW - eqn. 2.141

% Initialise matrices
Vn = cell(Nz,1);
inmGn = cell(Nz,1);
X_Upper = cell(Nz,1);
inmG_U = cell(Nz,1);
inmG_L = cell(Nz,1);
iexpG_U = cell(Nz,1);
expG_L = cell(Nz,1);
U = cell(Nz,1);
NullMat=zeros(2*Nt+1);
inmKX=diag(inmk0x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     The implementation of the stable algorithm begins
% (The notation is taken from Electromagnetic Surface Waves -
%  T. Mackay, J. Polo, and A. Lakhtakia, p75-77)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E=(Ex,Ey,Ez)=(Ex,0,Ez)+(0,Ey,0) and H=(Hx,Hy,Hz) = (Hx,0,Hz) +(0,Hy,0)
% Since E and H are orthogonal and by linearity, you can solve Maxwell
% for either E1=(Ex,0,Ez) and H1 = (0,Hy,0)
% or E2 = (0,Ey,0) and H2 = (Hx,0,Hz).
%(Ex,0,Ez)  and (0,Hy,0) corresponds to P-polaraization
%(0,Ey,0) and  (Hx,0,Hz) corresponds to S-polaraization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E = zeros(2*Nt+1, 2*Nt+1, Nz+1);
if sim.setup.Ez_A_or_iE
    iE = zeros(2*Nt+1, 2*Nt+1, Nz+1);
else
    A = zeros(2*Nt+1, 2*Nt+1, Nz+1);
end
%iA = zeros(2*Nt+1, 2*Nt+1, Nz+1);


epsz = 0;
iepsz = 0;

% Check flag for best method
if sim.setup.epsmethod == 2
    if sim.setup.pol == 1 % s-polarization use (1,0)
        epsmethod = 1;       
    else % p-polarization use (0,0)
        epsmethod = 0;       
    end
else
    epsmethod = sim.setup.epsmethod;
end


if sim.setup.iepsmethod == 2
    if sim.setup.pol == 1 % s-polarization use (1,0)
        iepsmethod = 0;       
    else % p-polarization use (0,0)
        iepsmethod = 0;       
    end
else
    iepsmethod = sim.setup.iepsmethod;
end


for n = Nz:-1:1
    % Select the optical permittivity, and the inverse, for the slice
    
    if ~isequal(epsz,epsf(n,:)) || ~isequal(iepsz,iepsf(n,:)) % See if slice has changed
        % update details
        epsz = epsf(n,:);
        iepsz = iepsf(n,:);
        
        % See if homogeneous
        neps = sum(abs(epsz)>=1e-12);
        nieps = sum(abs(iepsz)>=1e-12);
        %if neps == 1 || nieps == 1
        if 0
            Etemp =  epsz(2*Nt+1)*speye(2*Nt+1);
            Atemp = iepsz(2*Nt+1)*speye(2*Nt+1);
            iEtemp =  (1/epsz(2*Nt+1))*speye(2*Nt+1);
            iAtemp =  (1/iepsz(2*Nt+1))*speye(2*Nt+1);
            
            if sim.setup.Ez_A_or_iE
                iE(:,:,n) = iEtemp;
            else
                A(:,:,n) = Atemp;
            end
            
            if epsmethod
                eps = Etemp;
            else
                eps = iAtemp;
            end
            
            if iepsmethod
                ieps = Atemp;
            else
                ieps = iEtemp;
            end
                                    
            if (sim.setup.pol == 1)
                p_tr = sqrt(inmk0 - inmk0x.*inmk0x*(ieps(2*Nt+1,2*Nt+1)/inmk0)); % p[n]
                p_bl = sqrt(inmk0 * eps(2*Nt+1,2*Nt+1));    % p
                
                pd_p = diag(p_tr/p_bl);
                p_pd = diag(p_bl./p_tr);
            else
                p_tr = sqrt(-inmk0);
                p_bl = sqrt(- inmk0*eps(2*Nt+1,2*Nt+1) + (1./inmk0)*inmk0x.*inmk0x);
                
                pd_p = diag(p_tr./p_bl);
                p_pd = diag(p_bl/p_tr);
            end
            
            V =  [pd_p, -pd_p; eye(2*Nt+1), eye(2*Nt+1)];
            inmG = [p_bl*p_tr, -p_bl*p_tr];
            invV = 0.5*[p_pd eye(2*Nt+1) ; -p_pd eye(2*Nt+1)];
             
            [V, order] = sortmat(V,inmG);
            invV=invV(order,:);
            inmG = inmG(order);
            
        else
            % Create Toeplitz matrices - Eqn. 2.114
            Etemp = toeplitz(epsz(2*Nt+1:-1:1),epsz(2*Nt+1:4*Nt+1));  
            Atemp = toeplitz(iepsz(2*Nt+1:-1:1),iepsz(2*Nt+1:4*Nt+1));
            iEtemp = inv(Etemp);
            
            % This is needed later
            if sim.setup.Ez_A_or_iE
                iE(:,:,n) = inv(Etemp);
            else
                A(:,:,n) = Atemp;
            end
            
            % Take Toeplitz eps and ieps using different methods
            if epsmethod
                eps = Etemp;      
            else
                eps = inv(Atemp);
            end
            
            if iepsmethod
                ieps = Atemp;
            else
                if sim.setup.Ez_A_or_iE
                    ieps = iE(:,:,n);
                else
                    ieps = inv(Etemp);
                end
            end
            
            % Calculate P matrix components required for polarisation: 
            % ESW - Eqn. 2.119 with omege = inmk0 
            if sim.setup.pol == 1
                inmP14 = inmk0*eye(2*Nt + 1) - nmik0*inmKX*ieps*inmKX;
                inmP41 = inmk0*eps; 
                inmP = [NullMat, inmP14; inmP41, NullMat];
            elseif sim.setup.pol == 2
                inmP23 = - inmk0*eye(2*Nt + 1);
                inmP32 = - inmk0*eps + nmik0*inmKX*inmKX;           
                inmP = [NullMat, inmP23; inmP32, NullMat];
            else
                error('Polarization not defined');
            end
            
            % Perform eigendecomposition of P, such that P = V*G*iV 
            % Shown in ESW - eqn 2.138
            [V, inmGD] = eig(inmP);
            inmG = diag(inmGD);
                    

            % Sort the matrix by decreasing size of eigenvalue
            [V, order] = sortmat(V, inmG);
            inmG = inmG(order);
            
        end
    else

         % This is needed later
         if sim.setup.Ez_A_or_iE
                iE(:,:,n) = iEtemp;
            else
                A(:,:,n) = Atemp;
            end
            
    end
    
%     % Build X matrices ESW - Eqn. 2.144
%     if neps == 1 || nieps == 1
%         X = invV*Z{n+1};
%     else
%         X = V\Z{n+1}; % X = inv(V).Z{n+1} 
%     end
    Xc=rcond(V);
    if isnan(Xc)
        disp('NAN(V)')
        save junk.mat ll tt sim V inmP inmGD
    end
    if isnan(Xc) | Xc< 1.e-10
        X=pinv(V)*Z{n+1};
    else
        X=V\Z{n+1};
    end
    
    % store V and G in Vn and Nn for current layer
    Vn{n} = V;
    inmGn{n} = inmG;
    
    % Split G into upper and lower diagonal matrices
    inmG_U{n} = inmG(1:2*Nt + 1);
    inmG_L{n} = inmG(2*Nt + 2:4*Nt + 2);
    iexpG_U{n} = exp(1.i*nmdz(n).*inmG_U{n}); % Used in 2.146
    expG_L{n} = exp(-1.i*nmdz(n).*inmG_L{n}); % Used in 2.146
    
    % Break up the X matrix found in Eqn. 2.144
    X_Upper{n} = X(1:2*Nt+1, :);
    X_Lower{n} = X(2*Nt+2:4*Nt+2, :);
    
    % Define the U matrix as in Eqn. 2.146
    Xc=rcond(X_Upper{n});
    if isnan(Xc)
        disp('NAN(X_upper)')
    end
    if isnan(Xc) | Xc< 1.e-10
        U{n} = diag(expG_L{n}) * X_Lower{n} * (pinv(X_Upper{n}) * diag(iexpG_U{n}));
    else
        U{n} = diag(expG_L{n}) * X_Lower{n} * (X_Upper{n} \ diag(iexpG_U{n}));
    end

    % Calculate the next Z using Eqn. 2.145
    Z{n} = Vn{n}*[eye(2*Nt + 1); U{n}]; % Notice shift in Z as couting starts at 0
end

% We can now calculate Tn{i}
% A vector describing the excited modes of the incident field.
Ainc = zeros(2*Nt+1,1);  % Y and Z are for the correct pol only, so problem dims are halved
Ainc(Nt+1) = 1; % Zero mode excited

% Eqn. 2.149 typo in book, should read - [T0; R] = [Z(0), - Y_ref]^{-1} *
% Y_inc * Ainc
Xc=rcond([Z{1}, - Y_ref]);
if isnan(Xc)
    disp('T0R')
end
if isnan(Xc) | Xc < 1.e-10
    T0R =pinv([Z{1}, - Y_ref])*(Y_inc*Ainc);
else
    T0R = [Z{1}, - Y_ref]\(Y_inc*Ainc);
end
T0 = T0R(1 : 2*Nt + 1, :);
R = T0R(2*Nt + 2 : (4*Nt + 2), :);

% Build Transmission Vectors
Tn=cell(Nz+1,1);
Tn{1} = T0;

% Iteratively solve for Tn using 2.143 rearranged for T^(l)
for n = 1: Nz
    Xc=rcond(X_Upper{n});
    if isnan(Xc)
        disp('X_upper')
    end
    if isnan(Xc) | Xc < 1.e-10
        Tn{n+1} = pinv(X_Upper{n})*diag(iexpG_U{n})*Tn{n}; % Eqn. 2.143
    else
        Tn{n+1} = X_Upper{n}\diag(iexpG_U{n})*Tn{n}; % Eqn. 2.143
    end
end


% Calculate the electric fields
% Create the f vector which contains [e_x; e_y; h_x; h_y]
f = zeros(Nz,4*Nt+2);
f1 = zeros(Nz,4*Nt+2);
f2 = zeros(Nz,4*Nt+2);
e_x = zeros(Nz,2*Nt+1);
e_y = zeros(Nz,2*Nt+1);
e_z = zeros(Nz,2*Nt+1);
h_x = zeros(Nz,2*Nt+1);
h_y = zeros(Nz,2*Nt+1);
h_z = zeros(Nz,2*Nt+1);

for n = 1: Nz
    
     
     f1(n,:) = Z{n}*Tn{n};
     f2(n,:) = Z{n+1}*Tn{n+1};
     f(n,:) = 0.5*(f1(n,:)+f2(n,:)); % Eqn. 2.143
     
    if(sim.setup.pol == 1) % p-polarisation
        e_x(n,:) = f(n,1:2*Nt+1);
        h_y(n,:) = f(n,2*Nt+2:4*Nt+2);
        e_z(n,:) = - 1./inmk0 *  A(:,:,n) * (inmKX) * (h_y(n,:)).';       
     %  e_z(n,:) = - 1./inmk0 *  0.5 * (A(:,:,n) + A(:,:,n+1)) * (inmKX) * (h_y(n,:)).';       
     %  e_z(n,:) = -(0.5/inmk0)*(A(:,:,n) * (inmKX) * f1(n,2*Nt+2:4*Nt+2).'+ A(:,:,n+1) * (inmKX) * f2(n,2*Nt+2:4*Nt+2).');     
    elseif(sim.setup.pol == 2) % s-polarisation
        h_x(n,:) = f(n,2*Nt+2:4*Nt+2);
        e_y(n,:) = f(n,1:2*Nt+1);
        h_z(n,:) = (1./(inmk0 .* sim.phys.misc*sim.phys.mkgis2iA2mu0)) * inmKX*(e_y(n,:)).';
     end
end



% computations of transmission \
% coefficients, TP[i] contains the tranmission coefficients at the ith
% interface so TP[Nz+1] contains the coefficients of whole structure as
% their are Nz+1 interfaces
r2 = abs(R).^2;    % Reflection amplitudes
t2 = abs(Tn{end}).^2;    % Transmission amplitudes


% The following block computes a diagonal matrix with real parts of \
% kz^(n) as its elements**) ESW - 2.98 coefficients
kzr1 = real(inmk0z)/inmk0z(Nt+1);
RKZ1 = diag(kzr1);

%  (******************************* incident, reflected and \
% transmitted ***************************************************

R = r2;
T = RKZ1*t2;
R0 = R(Nt+1); % Specular reflection
R1 = sum(R)-R(Nt+1); % Scattered Reflection
T0 = T(Nt+1); % Specular Tranmission
T1 = sum(T)-T(Nt+1); % Scattered Transmission
% (***************************************************************


if sim.setup.pol == 2
    sim.zfl.Hx{ll} = sim.phys.Enormconst * h_x;
    sim.zfl.Hz{ll} = sim.phys.Enormconst * h_z;
    sim.zfl.Ey{ll} = sim.phys.Enormconst * e_y;
       
    sim.results.sR(ll) = abs(R0) + abs(R1);
    sim.results.sT(ll) = abs(T0) + abs(T1);
    sim.results.sA(ll) = 1 - sim.results.sR(ll) - sim.results.sT(ll);
    
elseif sim.setup.pol == 1
    
    sim.zfl.Ex{ll} = sim.phys.Enormconst * e_x;
    sim.zfl.Ez{ll} = sim.phys.Enormconst * e_z;
    sim.zfl.Hy{ll} = sim.phys.Enormconst * h_y;
   
    sim.results.pR(ll) = abs(R0) + abs(R1);
    sim.results.pT(ll) = abs(T0) + abs(T1);
    sim.results.pA(ll) = 1 - sim.results.pR(ll) - sim.results.pT(ll);
end

end




