function tests = ConversionsTest
% Unit tests for the leaf functions in Conversions/.
%
% Run locally with:
%   runtests('tests/ConversionsTest.m')
tests = functiontests(localfunctions);
end


% --- Length: nm <-> cm ---------------------------------------------------

function testCmFromNmKnownValue(testCase)
verifyEqual(testCase, cm_from_nm(1), 1e-7, 'RelTol', 1e-12);
end

function testNmFromCmKnownValue(testCase)
verifyEqual(testCase, nm_from_cm(1), 1e7, 'RelTol', 1e-12);
end

function testCmNmRoundtrip(testCase)
nm = [1, 250, 1064, 1e4];
verifyEqual(testCase, nm_from_cm(cm_from_nm(nm)), nm, 'RelTol', 1e-12);
end


% --- Length: nm <-> m ----------------------------------------------------

function testMFromNmKnownValue(testCase)
verifyEqual(testCase, m_from_nm(1), 1e-9, 'RelTol', 1e-12);
end

function testNmFromMKnownValue(testCase)
verifyEqual(testCase, nm_from_m(1), 1e9, 'RelTol', 1e-12);
end

function testMNmRoundtrip(testCase)
nm = [1, 250, 1064, 1e4];
verifyEqual(testCase, nm_from_m(m_from_nm(nm)), nm, 'RelTol', 1e-12);
end


% --- Inverse length: 1/nm -> 1/m -----------------------------------------

function testImFromInmKnownValue(testCase)
% 1 per nanometre = 1e9 per metre.
verifyEqual(testCase, im_from_inm(1), 1e9, 'RelTol', 1e-12);
end

function testImFromInmIsInverseOfMFromNm(testCase)
% (1/L) in inverse units should equal 1/(L in linear units).
nm = [1, 632.8, 1064];
verifyEqual(testCase, im_from_inm(1./nm), 1./m_from_nm(nm), 'RelTol', 1e-12);
end


% --- Angles --------------------------------------------------------------

function testRadFromDegKnownValues(testCase)
verifyEqual(testCase, rad_from_deg(0),   0,        'AbsTol', 1e-15);
verifyEqual(testCase, rad_from_deg(180), pi,       'RelTol', 1e-12);
verifyEqual(testCase, rad_from_deg(90),  pi/2,     'RelTol', 1e-12);
end

function testDegFromRadKnownValues(testCase)
verifyEqual(testCase, deg_from_rad(0),    0,   'AbsTol', 1e-12);
verifyEqual(testCase, deg_from_rad(pi),   180, 'RelTol', 1e-12);
verifyEqual(testCase, deg_from_rad(pi/2), 90,  'RelTol', 1e-12);
end

function testAngleRoundtrip(testCase)
deg = [-180, -45, 0, 30, 90, 180, 360];
verifyEqual(testCase, deg_from_rad(rad_from_deg(deg)), deg, 'RelTol', 1e-12);
end


% --- Photon energy from wavelength ---------------------------------------

function testEvFromNmKnownValues(testCase)
% hc = 1239.842 eV.nm in the module's chosen rounding.
verifyEqual(testCase, eV_from_nm(1239.842), 1.0, 'RelTol', 1e-9);

% Sanity: a 620 nm red photon is ~2 eV.
verifyEqual(testCase, eV_from_nm(620), 2.0, 'RelTol', 1e-3);

% Sanity: a 1240 nm IR photon is ~1 eV.
verifyEqual(testCase, eV_from_nm(1240), 1.0, 'RelTol', 1e-3);
end

function testEvFromNmIsVectorised(testCase)
nm = [400, 500, 620, 1240];
expected = 1239.842 ./ nm;
verifyEqual(testCase, eV_from_nm(nm), expected, 'RelTol', 1e-12);
end
