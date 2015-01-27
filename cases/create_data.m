function [y, u, A, B, c, lambda2_y, y_held] = create_data(N, M, D, noise)
% CREATE_DATA Creates a dataset and all test parameters for a CPCA problem
% using the init_true_values function, which sets A and B to have the same
% format as a D-tuple pendulum and samples from a von Mises using the 
% vmrand function. All samples are zero mean.

    % Generate u, A, B, c and lambda2_y from the init_true_val function
    [u, A, B, c, lambda2_y] = init_true_values(M, D, noise);
    C = repmat(c, [M/2, N]);
    
    % Sample from a von Mises
    theta = vmrand(0, repmat(u, [1, N]));

    % Generate the data from the model
    y = A * cos(theta) + B * sin(theta) + C + noise .* randn(M, N);
    
    N_held = fix(0.25 * N);
    theta_held = vmrand(0, repmat(u, [1, N_held]));
    C_held = repmat(c, [M/2, N_held]);
    y_held = A * cos(theta_held) + B * sin(theta_held) + ...
                C_held + noise .* randn(M, N_held);
end