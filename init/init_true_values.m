function [u, A, B, c, lambda2_y] = init_true_values(M, D, noise)
%INIT_TRUE_VALUES Initialisation at the true values of the model
    u = ones(D, 1) * 2.;
    A = zeros(M, D);
    B = zeros(M, D);
    % Create a structure similar to a D-tuple pendulum.
    for dd = 1:min(floor(M / 2), D)
        A(2 * dd - 1, 1:dd) = 1;
        B(2 * dd, 1:dd) = 1;
    end
    c = [0, 0]'; % no offset
    lambda2_y = noise^(-2);
end