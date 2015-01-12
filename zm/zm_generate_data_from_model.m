function [y, y_noiseless] = zm_generate_data_from_model(u, A, B, lambda2_y, d_pts)
%GENERATE_DATA_FROM_MODEL Generates data from the model
%   Detailed explanation goes here

    [m_dim, n_dim] = size(A);
    
    m_prior = zeros(size(u));
    k_prior = u;
    
    if n_dim == 1
        theta = vmrand(m_prior * ones(d_pts,1), k_prior * ones(d_pts,1))';
    else
        theta = vmrand(repmat(m_prior, [1, d_pts]), ...
                       repmat(k_prior, [1, d_pts]));
    end
    
    y_noiseless = A * cos(theta) + B * sin(theta);
    
    y = y_noiseless + 1 / lambda2_y * randn(m_dim, d_pts);
end

