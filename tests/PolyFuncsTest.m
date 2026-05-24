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


% --- lobpoly (Lagrange interpolants at LGL nodes) ------------------------

function testLobpolyKroneckerDeltaProperty(testCase)
% lobpoly(n) returns (n+1) Lagrange polynomials, where column j is the
% polynomial that equals 1 at the j-th LGL node and 0 at every other
% node. lobpoly flips the node order before building, so we evaluate
% against the same flipped grid.
n = 4;
P = lobpoly(n);
xi = fliplr(lgnodes(n));
for j = 1:n+1
    for k = 1:n+1
        if j == k
            verifyEqual(testCase, polyval(P(:,j), xi(k)), 1, ...
                'AbsTol', 1e-10);
        else
            verifyEqual(testCase, polyval(P(:,j), xi(k)), 0, ...
                'AbsTol', 1e-10);
        end
    end
end
end

function testLobpolyPartitionOfUnity(testCase)
% A nodal basis is a partition of unity: sum_j L_j(x) = 1 for any x.
n = 4;
P = lobpoly(n);
xs = linspace(-1, 1, 11);
for k = 1:numel(xs)
    total = 0;
    for j = 1:n+1
        total = total + polyval(P(:,j), xs(k));
    end
    verifyEqual(testCase, total, 1, 'AbsTol', 1e-10);
end
end


% --- node_av_global ------------------------------------------------------

function testNodeAvGlobalAveragesConstants(testCase)
% Two-element piecewise polynomial: constant 1 on the left, constant 2
% on the right. The single interior interface should be averaged to 1.5.
% Coefficient convention: column j is [a_n; ...; a_0] for element j.
left  = [0; 0; 1];   % polynomial == 1
right = [0; 0; 2];   % polynomial == 2
p = [left, right];
verifyEqual(testCase, node_av_global(p), 1.5, 'AbsTol', 1e-12);
end

function testNodeAvGlobalLinearOnEachElement(testCase)
% Three elements, each linear in xi in [-1, 1]:
%   p1(xi) = xi      -> right end value  1
%   p2(xi) = xi + 4  -> left end value   3, right end value  5
%   p3(xi) = 2*xi    -> left end value  -2
% Interface 1: 0.5*(1 + 3)   = 2
% Interface 2: 0.5*(5 + (-2)) = 1.5
p1 = [0; 1; 0];
p2 = [0; 1; 4];
p3 = [0; 2; 0];
p = [p1, p2, p3];
verifyEqual(testCase, node_av_global(p), [2; 1.5], 'AbsTol', 1e-12);
end
