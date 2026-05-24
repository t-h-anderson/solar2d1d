clear all
close all
clc;
addpath(genpath('./'))

% get default DE parameters
DEParams = DesignDE();

objFctSettings = 1;

% set title
optimInfo.title = 'Power Optimisation';

% specify objective function
objFctHandle = @OptoElec;


% define parameter names, ranges and quantization:

% 1. column: parameter names
% 2. column: parameter ranges
% 3. column: parameter quantizations
% 4. column: initial values (optional)

optical = 1;

switch 3
    case 1
        paramDefCell = {
            'nmLx', [100 1500], 1, 500;
            'nmLg', [1 500], 1, 150;
            'zeta', [0.01,0.99], 0.001, 0.5;
            
            'Eg0', [0.947, 1.626], 0.001, 0.95;
            'Eg1', [0.947, 1.626], 0.001, 1.5;
            
            'nmLi', [50, 1000], 1, 100;
            };
    case 2
        disp('case 2 under study')
        paramDefCell = {
% Tom suggests 
%             'nmLx', [100 1500], 1, 500;
%             'nmLg', [1 500], 1, 150;
%             'zeta', [0.01,0.99], 0.001, 0.5;
%             
%             'Eg_p', [0.947, 1.626], 0.001, 1;
%             'Eg_i', [0.947, 1.626], 0.001, 1;
%             'Eg_n', [0.947, 1.626], 0.001, 1;
%             
%             'nmLi', [50, 1000], 1, 100;
% Peter Main run
%            'nmLx', [100 1000], 20, 500;
%            'nmLg', [10 500], 10, 150;
%            'zeta', [0.01,0.99], 0.01, 0.5;          
%            'Eg_p', [0.95, 1.62], 0.01, 1.2;
%            'Eg_i', [0.95, 1.62], 0.01, 1.3;
%            'Eg_n', [0.95, 1.62], 0.01, 1.2;

% Peter refined after main           
            'nmLx', [250 800], 2, 620;
            'nmLg', [50 250], 1, 120;
            'zeta', [0.01,0.99], 0.01, 0.5;
            'Eg_p', [1.01, 1.62], 0.001, 1.4;      
            'nmLi', [100, 1000], 2 , 800;
            };
    case 3
                paramDefCell = {
            'nmLx', [100 1000], 1, 500;
            'nmLg', [1 500], 1, 150;
            'zeta', [0.01,0.99], 0.001, 0.5;
            
            'A', [0, 1], 0.0001, 1;
            'kappa', [0.01, 8], 0.01, 4;
            'phi', [0,1], 0.001, 0.75;
            'alpha', [-2,6], 0.001, 1;
            'Eg0', [0.947, 1.626], 0.001, 1;
            'nmLi', [50, 1000], 1, 100;
            };
end


% set times
DEParams.maxclock = [];

% set display options
DEParams.infoIterations = 1;
DEParams.infoPeriod = 2; % in seconds

% set random state in order to always use the same population members here
setrandomseed(1);

% start differential evolution
[bestmem, bestval, bestFctParams, nrOfIterations, resultFileName] = differentialevolution(DEParams, paramDefCell, objFctHandle, objFctSettings, [], [], optimInfo); %#ok

disp(' ');
disp('Best parameter set returned by function differentialevolution:');
disp(bestFctParams);

toexport=[optimResult.allTestedMembers;optimResult.allEvaluationValues;optimResult.allSubEvaluationValues];
csvwrite('./OptStudy.csv',toexport)
