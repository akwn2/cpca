function [u, v, A, B, lambda2_y] = init_true_val(m_dim, n_dim, useZeroMean, noise)
%INIT_TRUEVAL Initialisation by one pendulum at a time.
    
    A = zeros(m_dim, n_dim);
    B = zeros(m_dim, n_dim);
    
    for dd = 1:min(floor(m_dim / 2), n_dim)
        A(2 * dd - 1, 1:dd) = 1;
        B(2 * dd, 1:dd) = 1;
    end
    
    k_prior = ones(n_dim, 1) * 2.;
    m_prior = ones(n_dim, 1) * 0.;
    
    if useZeroMean
         u = k_prior;
         v = zeros(size(u));
    else
        u = k_prior .* cos(m_prior);
        v = k_prior .* sin(m_prior);
    end
    lambda2_y = noise^(-2);
end