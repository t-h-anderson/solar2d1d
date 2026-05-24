function tests = PolyFuncsTest
% Unit tests for the polynomial-helper leaves in DDNewton/PolyFuncs/.
%
% Coefficient convention follows MATLAB's polyval/polyder/conv: a column
% vector [a_n; a_{n-1}; ...; a_1; a_0] represents a_n*x^n + ... + a_0.
tests = functiontests(localfunctions);
end


% --- polytimes -----------------------------------------------------------

function testPolytimesNonzeroLeadingCoeffs(testCase)
% (x + 1)(x - 1) = x^2 - 1
p = [1; 1];
q = [1; -1];
expected = [1; 0; -1];
verifyEqual(testCase, polytimes(p, q), expected);
end

function testPolytimesShortcutsScalarP(testCase)
% polytimes spots when p degenerates to a constant and skips the convolution.
p = [0; 0; 5];   % constant 5 packed as a degree-2 polynomial
q = [1; 2; 3];   % x^2 + 2x + 3
verifyEqual(testCase, polytimes(p, q), 5 * q);
end

function testPolytimesShortcutsScalarQ(testCase)
% Symmetric case: q is a packed constant.
p = [1; 2; 3];
q = [0; 0; 7];
verifyEqual(testCase, polytimes(p, q), 7 * p);
end

function testPolytimesMatchesConv(testCase)
% For generic inputs polytimes should agree with conv.
rng(0);
p = randn(4, 1);
q = randn(3, 1);
verifyEqual(testCase, polytimes(p, q), conv(p, q), 'RelTol', 1e-12);
end


% --- polymatch -----------------------------------------------------------

function testPolymatchPadsShort(testCase)
% Pad the degree-1 polynomial 2x + 3 up to degree 3 by prepending zeros.
p = [2; 3];
verifyEqual(testCase, polymatch(p, 4), [0; 0; 2; 3]);
end

function testPolymatchTruncatesLong(testCase)
% Trim from the high-degree end (leading entries).
p = [1; 2; 3; 4; 5];
verifyEqual(testCase, polymatch(p, 3), [3; 4; 5]);
end

function testPolymatchIdentity(testCase)
% If the polynomial already has the requested length, leave it alone.
p = [4; 5; 6];
verifyEqual(testCase, polymatch(p, 3), p);
end


% --- pdiff ---------------------------------------------------------------

function testPdiffFirstDerivative(testCase)
% d/dx (x^3 + 2x + 1) = 3x^2 + 2; padded back to length 4.
p = [1; 0; 2; 1];
expected = [0; 3; 0; 2];
verifyEqual(testCase, pdiff(p, 1), expected);
end

function testPdiffSecondDerivative(testCase)
% d^2/dx^2 (x^3 + 2x + 1) = 6x; padded back to length 4.
p = [1; 0; 2; 1];
expected = [0; 0; 6; 0];
verifyEqual(testCase, pdiff(p, 2), expected);
end

function testPdiffOfConstantIsZero(testCase)
p = [0; 0; 7];
expected = [0; 0; 0];
verifyEqual(testCase, pdiff(p, 1), expected);
end


% --- lgnodes (Legendre-Gauss-Lobatto) -----------------------------------

function testLgnodesEndpointsAreClosed(testCase)
% LGL nodes always include the endpoints x = -1 and x = 1.
[x, w] = lgnodes(4);
verifyEqual(testCase, min(x), -1, 'AbsTol', 1e-12);
verifyEqual(testCase, max(x),  1, 'AbsTol', 1e-12);
% Weights for an LGL rule on [-1, 1] sum to the interval length, 2.
verifyEqual(testCase, sum(w), 2, 'RelTol', 1e-10);
end

function testLgnodesIntegratesPolynomialsExactly(testCase)
% N+1 LGL nodes integrate polynomials of degree up to 2N-1 exactly.
% For N=4 the integral of x^6 from -1 to 1 is 2/7.
[x, w] = lgnodes(4);
verifyEqual(testCase, sum(w .* x.^6), 2/7, 'RelTol', 1e-10);
% And the integral of x (odd) is 0.
verifyEqual(testCase, sum(w .* x), 0, 'AbsTol', 1e-12);
end
