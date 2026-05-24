function tests = LinAlgTest
% Unit tests for pure linear-algebra leaves in DDNewton/PolyFuncs/.
tests = functiontests(localfunctions);
end


% --- tiprod (rank-3 tensor x vector contraction) -------------------------

function testTiprodContractsDim1(testCase)
% Build a 2x3x4 tensor with X(i,j,k) = 100*i + 10*j + k. Contracting
% along dim 1 with Y = [a; b] yields a 3x4 matrix where entry (j,k) is
% a*(100 + 10j + k) + b*(200 + 10j + k).
X = zeros(2, 3, 4);
for i = 1:2
    for j = 1:3
        for k = 1:4
            X(i,j,k) = 100*i + 10*j + k;
        end
    end
end
Y = [1; 2];
expected = squeeze(X(1,:,:)) * Y(1) + squeeze(X(2,:,:)) * Y(2);
verifyEqual(testCase, tiprod(X, 1, Y), expected, 'RelTol', 1e-12);
end

function testTiprodContractsDim2(testCase)
% Contracting along dim 2 with a length-3 vector.
X = reshape(1:24, [2, 3, 4]);
Y = [1; 0; -1];
result = tiprod(X, 2, Y);
expected = squeeze(X(:,1,:)) * Y(1) + ...
           squeeze(X(:,2,:)) * Y(2) + ...
           squeeze(X(:,3,:)) * Y(3);
verifyEqual(testCase, result, expected, 'RelTol', 1e-12);
verifySize(testCase, result, [2, 4]);
end

function testTiprodContractsDim3(testCase)
% Contracting along dim 3 with a length-4 vector.
X = reshape(1:24, [2, 3, 4]);
Y = [1; 1; 1; 1];
result = tiprod(X, 3, Y);
expected = sum(X, 3);
verifyEqual(testCase, result, expected, 'RelTol', 1e-12);
verifySize(testCase, result, [2, 3]);
end

function testTiprodRejectsNonRank3Tensor(testCase)
% A rank-2 tensor must be rejected with a clear error.
X = ones(3, 3);
verifyError(testCase, @() tiprod(X, 1, [1; 1; 1]), ?MException);
end

function testTiprodRejectsMismatchedVectorLength(testCase)
% Y must have length equal to size(X, Xind).
X = ones(2, 3, 4);
verifyError(testCase, @() tiprod(X, 1, [1; 1; 1]), ?MException);
verifyError(testCase, @() tiprod(X, 2, [1; 1]),    ?MException);
verifyError(testCase, @() tiprod(X, 3, [1; 1; 1]), ?MException);
end

function testTiprodRejectsBadDimensionIndex(testCase)
X = ones(2, 3, 4);
verifyError(testCase, @() tiprod(X, 0, []), ?MException);
verifyError(testCase, @() tiprod(X, 4, []), ?MException);
end
