function sim = DesignSimOptDetail(sim)
% Design the simulation

sim.optical.setup = struct(...
     ... % Optical Parameters
    'optical_toggle', 0, ... % Toggle: 1 to run optical simulation, 0 to skip
    'parforArg', inf, ...   % Numer of parallel agents
    'epsmethod', 1, ...     % Formulaton used for eps in RCWA (2 For optimal, 1 for F(eps) and 0 (r) for invF(1/eps))
    'iepsmethod', 1, ...    % Formulaton used for ieps in RCWA (2 For optimal, 1 (r) for F(ieps) and 0 for invF(eps))
    'Ez_A_or_iE', 0, ...    % Choice for rebuilding Ez from Hy, 0 for A, 1 for iE
    'pol', 0, ...               % 0: Unpolarized, 1: p-polarized, 2: s-polarized
    ... % Analysis Options
    'outdim', 1, ... Dimension of output
    'avGdx_method' , 0, ... % Method of averaging G over x
    'RebuildFields', 0, ... % Use ifft to rebuild fields (overridden if Eplot or Hplot toggled)
    ... % Display options
    'optplot', 0, ...         % Produce plots (0: no, 1: yes)
    'epsplot', 1, ... Plot permittivity
    'Eplot', 1, ... Plot Efield
    'Hplot', 1, ... Plot Hfield
    ... % Mesh for angle and wavelength
    'degtheta0', 0, ...     % Minimum angle of Incident Light
    'degtheta1', 0, ...     % Maximum angle of Incident Light
    'ntheta', 1, ...        % Number of angles
    'nmlambda0', 680, ...   % Minimum wavelength in nm
    'nmlambda1', 680, ...  % Maximum wavelength in nm
    'nlambda', 1 + 0*floor((1100-350)/10), ...     % Number of wavelengths
    'Nx', 500, ...       % Number sample points in x-direction
    'nslices',  [200, 75, 20, 200, 20, 50, 100], ... % Number of slices in each sec (DesignJunction.m), uses 1 if not set.
    'dim', 2,...  % Output G in 1D (1) or 2D (2).     
    ... % Parameters for Nt convergence
    'CheckNtConvergence', 0, ... Use adaptive Nt
    'Nttol', 1e9,...      % Required tolerance in mA/(nm cm^2) for increasing Nt by 2
    'Nttol_lower', 1e-6, ...%  Required tolerance in mA/(nm cm^2) for decreasing Nt by 2 (used in serial calcs)
    'minNt', 20, ...         % Start value for Nt
    'maxNt', 20, ...        % Largest value for Nt
    'nosuccreq', 2, ... 
    'faillimit', 10, ...
    'ptol', 0.1 ...
);
          
sim.electrical.setup = struct(...
    ... % Electrical Parameters
    'electrical_toggle', 1, ... % Toggle: 1 to run electrical simulation, 0 to skip
    'pdeg', 3, ... % Polynomial degree
    'nx_sec', 20, ... % Mesh number
    'VdV', 0.05, ... % Change in volrage 
    ... % Simulation boundaries
    'mesh' , 1,...%-x.^2+x+scale*(1+scale) , ... % (Relative) lengths of each segment, scales to match Lz. Length either 1 or nx
    ..., % Global polys. The poly ax^2+bx+c -> [a, b, c], for example.
    'taun', 1, ...
    'taup', 1, ...
    'tauphi', 1, ...
    'upwind', 1, ...  % If toggled, will ignore value of t and upwind
    'tensor', 0, ... % Uses tensors internally
    ... % Local polys. If variable is set to inf then the global funct. (above) will be converted to the local function (below)
    'lp', 0, ...
    'tau_t', 1, ...
    't_1', 1e2, ...
    't_2', 1e-2, ...
    'quadtoggle', 0, ... Toggle quad vs lobatto integration
    'intdeg', 10, ... % Integration accuracy
    ... % Result plots
    'plotsolve', 0, ... % Number of steps to plot results. 0 for no output. Inf for first and last.
    'nx_plot', 10, ...
    'SolverSteps', 0, ... % Display solver steps
    'plot_densities', 1, ...
    'plot_fields', 1, ...
    'plot_levels', 1, ...
    'plot_currents', 1, ...
    'plot_U', 1, ...
    'plot_tau', 1, ...
  ... % Tolerance
    'dampmin', 0.1, ... % Set the initial damping level
    'damprelax', 1.2,... % Rate of relaxation of damping
    'trydamping', 3, ... % Number of times damping is increases to aid convergence
    'dampingincrease', 10, ... Amount to increase damping by each attempt
    'maxloop', 30, ... % Maximum number of iteration for convergence
    'reltol', 1e-2 ,... % Required relative tolerance
    'abstol', 1e-4, ... % Required absolute tolerance   
    'warn', 1, ... % Display warning messages
    'forcepositive', 0, ... % Force the carriers to be positive
    'test_J', 0 ... %Test the Jacobian
   );

%%%%%%%%%%%%%%%%%%%%%%
% Do not edit below this box %
%%%%%%%%%%%%%%%%%%%%%%

% Make sure nslices is correct lenght
diff =  sim.material.nsec - length(sim.optical.setup.nslices);
if diff > 0
    sim.optical.setup.nslices = padarray(sim.optical.setup.nslices, [0,1], diff, 'post');
end

sim.optical.setup.Nt = sim.optical.setup.minNt;

sim.electrical.setup.counter = 0;

% Create polynomial basis for simulation
sim.electrical.setup.poly = lobpoly(sim.electrical.setup.pdeg);
%sim.electrical.setup.poly=rand(sim.electrical.setup.pdeg1,sim.electrical.setup.pdeg1);
%sim.electrical.setup.poly = eye(sim.electrical.setup.pdeg1);

sim.electrical.setup.ipoly = inv(sim.electrical.setup.poly);


