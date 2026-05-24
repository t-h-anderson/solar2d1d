function tests = SimulationSnapshotTest
% Snapshot regression: a default-constructed solarcell.Simulation applied
% to a default Junction must materialise to a struct that is
% isequaln-equal to DesignSim(DesignJunction()).
%
% Two fields under electrical.setup are non-deterministic by construction:
%   start_time -- monotonic tic value, differs every call
%   total_time -- toc(start_time), tiny float, differs every call
% These are stripped from both sides before comparing.
tests = functiontests(localfunctions);
end


function testDefaultSimulationMatchesDesignSim(testCase)
expected = stripTimers(DesignSim(DesignJunction()));
actual   = stripTimers(solarcell.Simulation().toStruct(solarcell.Junction()));
verifyTrue(testCase, isequaln(actual, expected), ...
    formatSnapshotDiff(expected, actual));
end


function testOpticalSetupMatches(testCase)
expected = DesignSim(DesignJunction());
actual   = solarcell.Simulation().toStruct(solarcell.Junction());
verifyEqual(testCase, sort(fieldnames(actual.optical.setup)), ...
                       sort(fieldnames(expected.optical.setup)));
verifyTrue(testCase, isequaln(actual.optical.setup, expected.optical.setup), ...
    formatSectionDiff('optical.setup', expected.optical.setup, actual.optical.setup));
end


function testElectricalSetupMatchesModuloTimers(testCase)
expected = stripTimers(DesignSim(DesignJunction()));
actual   = stripTimers(solarcell.Simulation().toStruct(solarcell.Junction()));
verifyEqual(testCase, sort(fieldnames(actual.electrical.setup)), ...
                       sort(fieldnames(expected.electrical.setup)));
verifyTrue(testCase, isequaln(actual.electrical.setup, expected.electrical.setup), ...
    formatSectionDiff('electrical.setup', expected.electrical.setup, actual.electrical.setup));
end


function testJunctionFieldsCarryThrough(testCase)
% toStruct(junction) must leave the Junction-side fields untouched.
expected = DesignJunction();
actual   = solarcell.Simulation().toStruct(solarcell.Junction());
verifyEqual(testCase, actual.input, expected.input);
verifyTrue(testCase, isequaln(actual.material, expected.material));
verifyEqual(testCase, actual.optical.periodic, expected.optical.periodic);
verifyTrue(testCase, isequaln(actual.phys, expected.phys));
end


function testAcceptsRawJunctionStruct(testCase)
% Transitional: toStruct should also accept the legacy struct directly,
% not just a Junction object. Lets callers move over incrementally.
viaObject = stripTimers(solarcell.Simulation().toStruct(solarcell.Junction()));
viaStruct = stripTimers(solarcell.Simulation().toStruct(DesignJunction()));
verifyTrue(testCase, isequaln(viaObject, viaStruct));
end


function testNtInitialisedFromMinNt(testCase)
sim = solarcell.Simulation();
sim.minNt = 23;
s = sim.toStruct(solarcell.Junction());
verifyEqual(testCase, s.optical.setup.Nt, 23);
end


function testNslicesPaddedWhenShorterThanNsec(testCase)
% Verbatim fix-up from DesignSim.m: if nslices is shorter than nsec,
% it gets padded with ones. Pin that behaviour explicitly.
sim = solarcell.Simulation();
sim.nslices = [2, 3, 4];  % length 3 vs nsec=9
s = sim.toStruct(solarcell.Junction());
verifyEqual(testCase, s.optical.setup.nslices, [2, 3, 4, 1, 1, 1, 1, 1, 1]);
end


% --- helpers ------------------------------------------------------------

function s = stripTimers(s)
% Remove the two non-deterministic timer fields so snapshot comparison
% isn't defeated by sub-microsecond noise.
s.electrical.setup = rmfield(s.electrical.setup, {'start_time', 'total_time'});
end


function msg = formatSnapshotDiff(expected, actual)
fe = sort(fieldnames(expected));
fa = sort(fieldnames(actual));
lines = "solarcell.Simulation().toStruct(...) != DesignSim(DesignJunction()):";
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
