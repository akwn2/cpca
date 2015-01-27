function [u, A, B, c, lambda2_y] = init_perturb(M, D, noise, scale)
%INIT_PERTURB Initialisation by random perturbation of the true values with
% scaled random noise.

    % Get true values and perturbations
    [u_T, A_T, B_T, c_T, lambda2_y_T] = init_true_values(M, D, noise);
    [u_R, A_R, B_R, c_R, lambda2_y_R] = init_random(M, D);
    
    u = u_T + scale * u_R;
    A = A_T + scale * A_R;
    B = B_T + scale * B_R;
    c = c_T + scale * c_R;
    lambda2_y = lambda2_y_T + scale * lambda2_y_R;
end