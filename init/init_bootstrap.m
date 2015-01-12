function [u, v, A, B, lambda2_y] = init_bootstrap(y, u, v, A, B, lambda2_y, ...
    y_held, useZeroMean)
%DO_INIT_NTH_PENDULUM Initialisation by one pendulum at a time.

if nargin < 7
    useZeroMean = 0;
end
%   Detailed explanation goes here
    n_dim = size(A, 2);
    [m_dim, d_pts] = size(y);
    
    % FIXME - CHANGE TO GET APPROPRIATE DIMENSIONS
    for ii = 1:min(n_dim, m_dim / 2)
        k_prior = abs(u(ii) + 1.i * v(ii));
        m_prior = angle(u(ii) + 1.i * v(ii));

        k1 = k_prior * ones(1, d_pts);
        k2 = zeros(1, d_pts);
        m1 = m_prior * ones(1, d_pts);
        m2 = zeros(1, d_pts);

        [new_u, new_v, new_A, new_B, ~, ~] = ...
            do_vem(y(2 * ii-1:2 * ii, :), u(ii), v(ii), ...
            A(2*ii-1:2*ii, ii), B(2*ii-1:2*ii, ii), lambda2_y, ...
            k1, k2, m1, m2);

        u(ii) = new_u;
        v(ii) = new_v;
        A(2 * ii-1:2 * ii, ii) = new_A;
        B(2 * ii-1:2 * ii, ii) = new_B;
%         keyboard;
    end
end

