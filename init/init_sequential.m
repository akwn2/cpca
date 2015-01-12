function [u, v, A, B, lambda2_y] = init_sequential(y, u, v, A, B, lambda2_y)
%INIT_SEQUENTIAL Initialisation by building a single pendulum on
%top of learned pendulums.
%   Detailed explanation goes here
    n_dim = size(A, 2);
    [m_dim, d_pts] = size(y);
    
    for ii = 1:min(n_dim, m_dim / 2)
        
        % Constrain to submatrix
        Asub = A(1:2 * ii, 1:ii);
        Bsub = B(1:2 * ii, 1:ii);
        ysub = y(1:2 * ii, :);
        
        usub = u(1:ii);
        vsub = v(1:ii);
        
        k_prior = abs(usub + 1.i * vsub);
        m_prior = angle(usub + 1.i * vsub);

        k1 = repmat(k_prior, [1, d_pts]);
        k2 = zeros(size(k1));
        m1 = repmat(m_prior, [1, d_pts]);
        m2 = zeros(size(m1));

        % Learn
        [new_u, new_v, new_A, new_B, ~, ~] = ...
            init_sequential_do_vem(ysub, usub, vsub, ...
                                   Asub, Bsub, lambda2_y, ...
                                   k1, k2, m1, m2);

        % Update from the submatrix
        u(1:ii) = new_u;
        v(1:ii) = new_v;
        A(1:2 * ii, 1:ii) = new_A;
        B(1:2 * ii, 1:ii) = new_B;
%         keyboard;
    end
end

