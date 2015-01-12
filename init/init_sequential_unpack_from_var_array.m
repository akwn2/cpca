function [u, v, A, B, lambda2_y] = init_sequential_unpack_from_var_array(var_array, m_dim)
%UNPACK_FROM_VAR_ARRAY Unpacks model parameters from variable array
%   Detailed explanation goes here
    u = var_array(1);
    v = var_array(2);
    A = var_array(3: 2 + m_dim);
    B = var_array(3 + m_dim: 2 + 2 * m_dim);
    lambda2_y = var_array(3 + 2 * m_dim);
    
    A = reshape(A, [m_dim, 1]);
    B = reshape(B, [m_dim, 1]);
end

