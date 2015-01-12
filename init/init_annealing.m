function [u, v, A, B, lambda2_y] = init_annealing(y, u, v, A, B, ...
    y_held, useZeroMean, coolingPasses, coolingFactor)
%INIT_ANNEALING Initialisation by annealing

if nargin < 7
    useZeroMean = 0;
end
if nargin < 8
    coolingPasses = 5;
end
if nargin < 9
    coolingFactor = 2;
end

d_pts = size(y, 2);
n_dim = size(A, 2);

if useZeroMean
    k1 = repmat(u, [1, d_pts]);
    k2 = small_rand(n_dim, d_pts);
    m1 = small_rand(n_dim, d_pts);
    m2 = small_rand(n_dim, d_pts);
else
    k_prior = abs(u + 1.i * v);
    m_prior = angle(u + 1.i * v);
    
    k1 = repmat(k_prior, [1, d_pts]);
    k2 = small_rand(n_dim, d_pts);
    m1 = repmat(m_prior, [1, d_pts]);
    m2 = small_rand(n_dim, d_pts);
end

sigma2 = 1000;
lambda2_y = 1/sigma2;

for ll = 1:coolingPasses
    
    if useZeroMean
            lambda2_y = lambda2_y / ll;
            [u, A, B, lambda2_y, ~] = ...
                fixL_do_vem(y, u, A, B, lambda2_y, k1, k2, m1, m2, y_held);
    else
        [u, v, A, B, lambda2_y, ~] = ...
            do_vem(y, u, v, A, B, lambda2_y, k1, k2, m1, m2, y_held);
    end
    
    lambda2_y = lambda2_y / coolingFactor;
end