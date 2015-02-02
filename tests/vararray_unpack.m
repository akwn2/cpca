function [u, A, B, lambda2_y] = vararray_unpack(vararray, D, M)
%UNPACK_FROM_VAR_ARRAY Unpacks model parameters from variable array
%   Detailed explanation goes here
    u = vararray(1:D);
    A = vararray(D + 1: D + D * M);
    B = vararray(D + D * M + 1: D + 2 * (D * M));
    lambda2_y = vararray(D + 2 * (D * M) + 1);
    
    A = reshape(A, [M, D]);
    B = reshape(B, [M, D]);
end

