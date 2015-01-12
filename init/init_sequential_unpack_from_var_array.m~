function [u, v, A, B, lambda2_y] = unpack_from_var_array(var_array, n_dim, m_dim)
%UNPACK_FROM_VAR_ARRAY Unpacks model parameters from variable array
%   Detailed explanation goes here
    u = var_array(1:n_dim);
    v = var_array(n_dim + 1:2*n_dim);
    A = var_array(2 * n_dim + 1: 2 * n_dim + n_dim * m_dim);
    B = var_array(2 * n_dim + n_dim * m_dim + 1: 2 * (n_dim + n_dim * m_dim));
    lambda2_y = var_array(2 * (n_dim + n_dim * m_dim) + 1);
    
    A = reshape(A, [m_dim, n_dim]);
    B = reshape(B, [m_dim, n_dim]);
end

