function tests = CodeStyleTest
% Structural / hygiene checks that run alongside the regular unit tests.
tests = functiontests(localfunctions);
end


function testFunctionFilesDeclareMatchingFunctionName(testCase)
% Every .m file that opens with a `function` declaration in one of the
% active source directories should declare a function whose name matches
% the filename. Script files (no leading `function`) are skipped.
%
% Catches copy-paste residues like Conversions/cm_from_nm.m once
% declaring `function ... = cm2nm(...)`, or Conversions/ibp2sbp.m
% declaring its inner function as `sbp2ibp`. MATLAB resolves by
% filename so these compile fine, but they lie about what they compute
% and break reflection-based tooling.

repoRoot = fileparts(fileparts(mfilename('fullpath')));

% Active source directories. legacy/, tests/, DEA/ (third-party),
% Slaves/ (transient), Saved Inputs/ (configuration), Materials data
% subfolders, and Mathematica notebooks are all out of scope.
scanDirs = { ...
    '.', ...
    'Conversions', ...
    'DDNewton', ...
    fullfile('DDNewton', 'Display'), ...
    fullfile('DDNewton', 'PolyFuncs'), ...
    'RCWA', ...
    'Materials'};

mismatches = strings(0, 1);
filesScanned = 0;
functionFilesChecked = 0;

for d = 1:numel(scanDirs)
    folder = fullfile(repoRoot, scanDirs{d});
    if exist(folder, 'dir') ~= 7
        continue;  % directory removed since the test was written
    end
    files = dir(fullfile(folder, '*.m'));
    for k = 1:numel(files)
        filesScanned = filesScanned + 1;
        [~, expected] = fileparts(files(k).name);
        info = readFirstFunctionName(fullfile(files(k).folder, files(k).name));
        if ~info.isFunctionFile
            continue;  % plain script -- nothing to check
        end
        functionFilesChecked = functionFilesChecked + 1;
        if ~strcmp(info.name, expected)
            mismatches(end+1, 1) = sprintf('%s/%s declares "%s"', ...
                scanDirs{d}, files(k).name, info.name); %#ok<AGROW>
        end
    end
end

verifyGreaterThan(testCase, filesScanned, 0, ...
    'Guard test scanned zero files -- scanDirs likely misconfigured');
verifyGreaterThan(testCase, functionFilesChecked, 0, ...
    'Guard test found no function files -- parser likely broken');

verifyEmpty(testCase, mismatches, sprintf( ...
    'Function name does not match filename for:\n  %s', ...
    strjoin(mismatches, sprintf('\n  '))));
end


function info = readFirstFunctionName(filePath)
% Return a struct describing the first declaration in filePath.
%   info.isFunctionFile  true iff the first non-comment, non-blank line
%                        starts with `function`
%   info.name            the declared function name, or '' if the file
%                        is a script or the declaration could not be
%                        parsed
info = struct('isFunctionFile', false, 'name', '');

fid = fopen(filePath, 'r');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

while true
    line = fgetl(fid);
    if ~ischar(line)
        return;
    end
    stripped = strtrim(line);
    if isempty(stripped) || stripped(1) == '%'
        continue;
    end
    if isempty(regexp(stripped, '^function\s', 'once'))
        return;  % script file (e.g. starts with `clear`, `addpath`, ...)
    end
    info.isFunctionFile = true;
    % Match: function [outs] = name(args)  |  function name(args)  |
    %        function out = name           |  function name
    tokens = regexp(stripped, '^function\s+(?:[^=]*=\s*)?([A-Za-z_]\w*)', ...
        'tokens', 'once');
    if ~isempty(tokens)
        info.name = tokens{1};
    end
    return;
end
end
