function tests = CodeStyleTest
% Structural / hygiene checks that run alongside the regular unit tests.
tests = functiontests(localfunctions);
end


function testConversionsFunctionNamesMatchFilenames(testCase)
% Guard against copy-paste residues like the historical cm_from_nm.m and
% nm_from_cm.m, both of which declared `function ... = cm2nm(...)` even
% though one performed the reverse conversion. MATLAB resolves by
% filename so this compiles fine, but it breaks any reflection-based
% tooling and lies to the reader.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
folder = fullfile(repoRoot, 'Conversions');
files = dir(fullfile(folder, '*.m'));

verifyNotEmpty(testCase, files, ...
    'Conversions/ contains no .m files -- guard is meaningless');

mismatches = strings(0, 1);
for k = 1:numel(files)
    [~, expected] = fileparts(files(k).name);
    declared = readFirstFunctionName(fullfile(files(k).folder, files(k).name));
    if isempty(declared)
        mismatches(end+1, 1) = sprintf('%s: no function declaration parsed', ...
            files(k).name); %#ok<AGROW>
    elseif ~strcmp(declared, expected)
        mismatches(end+1, 1) = sprintf('%s declares "%s"', ...
            files(k).name, declared); %#ok<AGROW>
    end
end

verifyEmpty(testCase, mismatches, sprintf( ...
    'Function name does not match filename for:\n  %s', ...
    strjoin(mismatches, sprintf('\n  '))));
end


function name = readFirstFunctionName(filePath)
% Return the name of the first non-comment function declaration in
% filePath, or '' if none is found.
name = '';
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
    % Match: function [outs] = name(args)  |  function name(args)  |
    %        function out = name           |  function name
    tokens = regexp(stripped, '^function\s+(?:[^=]*=\s*)?([A-Za-z_]\w*)', ...
        'tokens', 'once');
    if ~isempty(tokens)
        name = tokens{1};
    end
    return;
end
end
