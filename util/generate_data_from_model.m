function y = generate_data_from_model(u, v, A, B, lambda2_y, d_pts)
%GENERATE_DATA_FROM_MODEL Generates data from the model
%   Detailed explanation goes here

    [m_dim, n_dim] = size(A);
    
    m_prior = angle(u + 1.i * v);
    k_prior = abs(u + 1.i * v);
    
    if n_dim == 1
        theta = vmrand(m_prior * ones(d_pts,1), k_prior * ones(d_pts,1))';
    else
        theta = vmrand(repmat(m_prior, [1, d_pts]), ...
                       repmat(k_prior, [1, d_pts]));
    end
    
    y = A * cos(theta) + B * sin(theta);
    
%     y = y_noiseless + 1 / lambda2_y * randn(m_dim, d_pts);
end

