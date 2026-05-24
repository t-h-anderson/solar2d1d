function DEParamsDefault = DesignDE
%GETDEFAULTPARAMS  Get default parameters for differential evolution.
%		DEParams = GETDEFAULTPARAMS returns a structure with a set of default
%		parameters for differential evolution.
%
%		<a href="differentialevolution.html">differentialevolution.html</a>  <a href="http://www.mathworks.com/matlabcentral/fileexchange/18593">File Exchange</a>  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=KAECWD2H7EJFN">Donate via PayPal</a>
%
%		Markus Buehren
%		Last modified 05.07.2011
%
%		See also DIFFERENTIALEVOLUTION.

% Basic DE parameters
DEParamsDefault.algorithm      = 'DE'; % optimization algorithm
DEParamsDefault.VTR            = NaN;  % "value to reach"
DEParamsDefault.NP             = 100;   % number of population members
DEParamsDefault.F              = -1;    % DE-stepsize <=0 for dither
DEParamsDefault.CR             = 0.8;  % crossover probability
DEParamsDefault.strategy       = 2;    % DE strategy

% parameters for finishing the optimization
DEParamsDefault.maxiter        = inf;  % maximum number of iterations
DEParamsDefault.maxtime        = inf;  % maximum time in seconds
DEParamsDefault.maxclock       = [];   % time when to finish optimization (time vector as 
%DEParamsDefault.maxtime        = 365*7*24*60*60; % 1 year in seconds                                       % returned by clock.m)

% parameters for reinitialization of population
DEParamsDefault.minvalstddev   = -1;
DEParamsDefault.minparamstddev = -1;
DEParamsDefault.nofevaliter    = 10;
DEParamsDefault.nochangeiter   = 10;

% parameters for information display
DEParamsDefault.infoIterations = 10;   % number of iterations
DEParamsDefault.infoPeriod     = 60;   % in seconds
DEParamsDefault.sendMailPeriod = 1800; % in seconds

% slave process parameters
DEParamsDefault.feedSlaveProc  = 1;
DEParamsDefault.slaveFileDir   = strcat(pwd, '/Slaves');
DEParamsDefault.maxMasterEvals = inf;

% miscellaneous
DEParamsDefault.useInitParams  = 1;
DEParamsDefault.saveHistory    = 1;
DEParamsDefault.displayResults = 0; % Leave zero - crashes otherwise
DEParamsDefault.playSound      = 1;
DEParamsDefault.minimizeValue  = 0;
DEParamsDefault.validChkHandle = '';

% Set up the list of returned variables from the objective funtion. Must
% manually add these in function, e.g. in OptoElec.m
DEParamsDefault.subFoMNames={'JscOpt', 'JscOpts', 'JscOptp', 'eta', 'Vmax', 'Jsc', 'Voc', 'FF','Optical Time', 'Electrical Time'};
DEParamsDefault.subFoM = length(DEParamsDefault.subFoMNames);



