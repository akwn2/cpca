function [u, v, A, B, lambda2_y] = init_randomly(m_dim, n_dim, useZeroMean, w, large_s2)

    if nargin < 3
        useZeroMean = 0;
    end
    if nargin < 4
        w = 1E0;
    end
    if nargin < 5
        large_s2 = 1E2;
    end

    A = rand(m_dim, n_dim) / w;
    B = rand(m_dim, n_dim) / w;

    k_prior = rand(n_dim, 1) / w;
    m_prior = zeros(n_dim, 1) / w;

    if useZeroMean
        u = k_prior;
        v = zeros(size(u));
    else
        u = k_prior .* cos(m_prior);
        v = k_prior .* sin(m_prior);
    end

    lambda2_y = 1 / large_s2;
end