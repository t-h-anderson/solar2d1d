function results = runAllTests(varargin)
% Discover and run every matlab.unittest case under tests/.
%
% Adds the repo root (and everything below it) to the path so test files
% can call into Conversions/, DDNewton/, RCWA/, etc. without the caller
% having to set things up first.
%
% Usage:
%   results = runAllTests                       % local interactive
%   results = runAllTests('JUnit', 'out.xml')   % write JUnit XML
%   results = runAllTests('Cobertura', 'cov.xml') % write coverage XML
%   results = runAllTests('FailOnError', true)  % error() on any failure
%
% Returns the matlab.unittest.TestResult array.

p = inputParser;
p.addParameter('JUnit', '', @ischar);
p.addParameter('Cobertura', '', @ischar);
p.addParameter('FailOnError', false, @islogical);
p.parse(varargin{:});
opts = p.Results;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

suite = matlab.unittest.TestSuite.fromFolder( ...
    fullfile(repoRoot, 'tests'), 'IncludingSubfolders', true);

runner = matlab.unittest.TestRunner.withTextOutput;

if ~isempty(opts.JUnit)
    ensureFolder(opts.JUnit);
    runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat(opts.JUnit));
end

if ~isempty(opts.Cobertura)
    ensureFolder(opts.Cobertura);
    sourceFolders = { ...
        fullfile(repoRoot, 'Conversions'), ...
        fullfile(repoRoot, 'DDNewton'), ...
        fullfile(repoRoot, 'RCWA'), ...
        fullfile(repoRoot, 'Materials')};
    sourceFolders = sourceFolders(cellfun(@(d) exist(d, 'dir') == 7, sourceFolders));
    plugin = matlab.unittest.plugins.CodeCoveragePlugin.forFolder(sourceFolders, ...
        'IncludingSubfolders', true, ...
        'Producing', matlab.unittest.plugins.codecoverage.CoberturaFormat(opts.Cobertura));
    runner.addPlugin(plugin);
end

results = runner.run(suite);

if opts.FailOnError && (any([results.Failed]) || any([results.Incomplete]))
    error('runAllTests:Failures', '%d failed, %d incomplete', ...
        nnz([results.Failed]), nnz([results.Incomplete]));
end

end


function ensureFolder(filePath)
folder = fileparts(filePath);
if ~isempty(folder) && exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
