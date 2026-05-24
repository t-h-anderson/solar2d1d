function tests = JunctionSnapshotTest
% Snapshot regression: a default-constructed solarcell.Junction must
% materialise to a struct that is isequaln-equal to DesignJunction()'s
% return value. This pins behaviour while we peel callers off the legacy
% function -- if the class drifts, this test catches it.
tests = functiontests(localfunctions);
end


function testDefaultJunctionMatchesDesignJunction(testCase)
expected = DesignJunction();
actual   = solarcell.Junction().toStruct();
verifyTrue(testCase, isequaln(actual, expected), ...
    formatSnapshotDiff(expected, actual));
end


function testInputSectionMatches(testCase)
% Per-section asserts give a more useful failure message when the
% top-level snapshot diverges.
expected = DesignJunction();
actual   = solarcell.Junction().toStruct();
verifyEqual(testCase, fieldnames(actual.input), fieldnames(expected.input));
verifyEqual(testCase, actual.input, expected.input);
end


function testMaterialSectionMatches(testCase)
expected = DesignJunction();
actual   = solarcell.Junction().toStruct();
verifyEqual(testCase, sort(fieldnames(actual.material)), ...
                       sort(fieldnames(expected.material)));
verifyTrue(testCase, isequaln(actual.material, expected.material));
end


function testOpticalPeriodicMatches(testCase)
expected = DesignJunction();
actual   = solarcell.Junction().toStruct();
verifyEqual(testCase, actual.optical.periodic, expected.optical.periodic);
end


function testPhysSectionMatches(testCase)
% phys is the trickiest section because most fields are derived from
% material+input data. isequaln so any NaN propagation in either side
% compares equal.
expected = DesignJunction();
actual   = solarcell.Junction().toStruct();
verifyEqual(testCase, sort(fieldnames(actual.phys)), ...
                       sort(fieldnames(expected.phys)));
verifyTrue(testCase, isequaln(actual.phys, expected.phys), ...
    formatSectionDiff('phys', expected.phys, actual.phys));
end


function testValidationRejectsBadStackLength(testCase)
% nsec=9 but nmSec only length 3 -> validate() must error.
j = solarcell.Junction();
j.nsec = 9;
j.nmSec = [100, 200, 300];
verifyError(testCase, @() j.toStruct(), ...
    'solarcell:Junction:lengthMismatch');
end


% --- diagnostics --------------------------------------------------------

function msg = formatSnapshotDiff(expected, actual)
% Build a human-readable summary of which top-level fields disagree.
fe = sort(fieldnames(expected));
fa = sort(fieldnames(actual));
lines = "solarcell.Junction().toStruct() != DesignJunction():";
missing = setdiff(fe, fa);
extra   = setdiff(fa, fe);
if ~isempty(missing)
    lines(end+1, 1) = sprintf('  missing in actual: %s', strjoin(missing, ', '));
end
if ~isempty(extra)
    lines(end+1, 1) = sprintf('  unexpected in actual: %s', strjoin(extra, ', '));
end
for k = 1:numel(fe)
    f = fe{k};
    if isfield(actual, f) && ~isequaln(actual.(f), expected.(f))
        lines(end+1, 1) = sprintf('  field "%s" differs', f); %#ok<AGROW>
    end
end
msg = char(strjoin(lines, sprintf('\n')));
end


function msg = formatSectionDiff(name, expected, actual)
fe = sort(fieldnames(expected));
fa = sort(fieldnames(actual));
lines = string(sprintf('Section "%s" differs:', name));
missing = setdiff(fe, fa);
extra   = setdiff(fa, fe);
if ~isempty(missing)
    lines(end+1, 1) = sprintf('  missing: %s', strjoin(missing, ', ')); %#ok<AGROW>
end
if ~isempty(extra)
    lines(end+1, 1) = sprintf('  unexpected: %s', strjoin(extra, ', ')); %#ok<AGROW>
end
for k = 1:numel(fe)
    f = fe{k};
    if isfield(actual, f) && ~isequaln(actual.(f), expected.(f))
        lines(end+1, 1) = sprintf('  "%s": expected %s, got %s', ...
            f, mat2str(expected.(f)), mat2str(actual.(f))); %#ok<AGROW>
    end
end
msg = char(strjoin(lines, sprintf('\n')));
end
