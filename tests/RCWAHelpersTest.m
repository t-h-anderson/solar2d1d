function tests = RCWAHelpersTest
% Unit tests for the stand-alone helpers in RCWA/.
tests = functiontests(localfunctions);
end


% --- sortmat -------------------------------------------------------------

function testSortmatOrdersColumnsByImagOfDuDescending(testCase)
% sortmat permutes the columns of G in the same order that sorting
% imag(Du) in descending order produces.
G = [10, 20, 30, 40;
     11, 22, 33, 44];
Du = [-1; 2i; 1+5i; 0.5];      % imag = [0, 2, 5, 0]
expectedOrder = [3; 2; 1; 4];  % imag sorted desc: 5, 2, 0, 0
[GSorted, Isort] = sortmat(G, Du);
verifyEqual(testCase, Isort, expectedOrder);
verifyEqual(testCase, GSorted, G(:, expectedOrder));
end

function testSortmatPreservesShapeOfRealMatrix(testCase)
% With purely real eigenvalues the order is determined by the (tied) zero
% imaginary part, but the output must still have the same shape as G.
G = magic(4);
Du = [1; 2; 3; 4];
[GSorted, Isort] = sortmat(G, Du);
verifySize(testCase, GSorted, size(G));
verifySize(testCase, Isort, [4, 1]);
% The permutation must be a permutation of 1:4.
verifyEqual(testCase, sort(Isort), (1:4).');
end
