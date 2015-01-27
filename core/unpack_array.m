function [u, A, B, c, lambda2_y] = unpack_array(var_array, D, M)
%UNPACK_ARRAY Unpacks model parameters from variable array
%   Detailed explanation goes here
    u = var_array(1:D);
    A = var_array(D + 1: D + D * M);
    B = var_array(D + D * M + 1: D + 2 * (D * M));
    c = var_array(D + 2 * (D * M) + 1: D + 2 * (D * M) + 2);
    lambda2_y = var_array(D + 2 * (D * M) + 3);
    
    A = reshape(A, [M, D]);
    B = reshape(B, [M, D]);
end

