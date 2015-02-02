function [u, A, B, c, lambda2_y, k1, k2, m1, m2] = opt_free_energy(y,...
                                    u, A, B, c, lambda2_y, ...
                                    k1, k2, m1, m2)
%OPT_FREE_ENERGY Learn an CPCA model by directly optimising the model free
% energy.
   
    % Get problem dimensions
    [M, D] = size(A);
    N = size(k1, 2);
    
    % Reshape and pack for optimization
    A = reshape(A, [M * D, 1]);
    B = reshape(B, [M * D, 1]);
    
    k1 = reshape(k1, [N * D, 1]);
    k2 = reshape(k2, [N * D, 1]);
    m1 = reshape(m1, [N * D, 1]);
    m2 = reshape(m2, [N * D, 1]);
    
    vararray = cat(1, u, A, B, c, lambda2_y, k1, k2, m1, m2);
    
    % Optimizing
    vararray = minimize(vararray, 'negative_free_energy', - max_iter, y, N, M, D);
    
    % Reshaping to output
    
    % Model parameters
    u = vararray(1:D);
    A = reshape(vararray(D + 1: D + D * M), [M, D]);
    B = reshape(vararray(D + D * M + 1: D + 2 * D * M), [M, D]);
    c = vararray(D + 2 * D * M + 1: D + 2 * D * M + 2);
    lambda2_y = vararray(D + 2 * D * M + 3);
    
    % Mean field terms
    k1 = reshape(vararray(D + 2 * D * M + 4:D + 2 * D * M + D * N + 3), [D, N]);
    k2 = reshape(vararray(D + 2 * D * M + D * N + 4:D + 2 * D * M + 2 * D * N + 3), [D, N]);
    m1 = reshape(vararray(D + 2 * D * M + 2 * D * N + 4 : D + 2 * D * M + 3 * D * N + 3), [D, N]);
    m2 = reshape(vararray(D + 2 * D * M + 3 * D * N + 4 : D + 2 * D * M + 4 * D * N + 3), [D, N]);
    
end