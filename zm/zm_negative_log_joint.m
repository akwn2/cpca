function [log_p, dlog_p] = negative_log_joint(var_array, y, mc, ms, mc_sq, ms_sq, msc, n_dim, m_dim, d_pts)
%NEGATIVE_LOG_JOINT Negative of the log-joint of for the M-step
%   Negative log joint and gradients for the probabilistic model of CPCA
    
    [log_u, A, B, log_lambda2_y] = zm_unpack_from_var_array(var_array, n_dim, m_dim);
    
    % General pre-computable relations
    lambda2_y = exp(log_lambda2_y);
    u = exp(log_u);
    
    U = repmat(u, [1, d_pts]);
    YYT = y * y';
    
    ATA = A' * A;
    ATB = A' * B;
    BTB = B' * B;
    
    sum_ccT = mc * mc' + diag(- diag(mc * mc') + sum(mc_sq, 2));
    sum_scT = ms * mc' + diag(- diag(ms * mc') + sum(msc, 2));
    sum_ssT = ms * ms' + diag(- diag(ms * ms') + sum(ms_sq, 2));
    
    % Gradient of the average log joint calculation
    log_p = d_pts * (n_dim - m_dim / 2.) * log(2. * pi) ... %               2 pi constant from prior and likelihood
        - d_pts * sum(u + log(besseli(0, u, 1))) ... %                      Prior normalizing constant
        + trace(U' * mc) ... %                                              Prior "energy"
        + 0.5 * d_pts * m_dim * log(lambda2_y) ... %                         Likelihood normalizing constant
        - 0.5 * lambda2_y * (trace(YYT) ... %                                   sum_{m, d} y_{m,d} y_{m,d}
                                - 2. * trace(mc' * A' * y) ... %                sum_{n, m, d} (mc_{d,n} A_{n,m}) y_{m,d}
                                - 2. * trace(ms' * B' * y) ... %                sum_{n, m, d} (ms_{d,n} B_{n,m}) y_{m,d}
                                + trace(ATA * sum_ccT) ... %                    sum_{n, m, d} (A_{n,m} A_{m,n}) (mc_{n,d} mc_{d,n})
                                + 2 * trace(ATB * sum_scT) ... %                sum_{n, m, d} (A_{n,m} B_{m,n}) (ms_{n,d} mc_{d,n})
                                + trace(BTB * sum_ssT)); %                      sum_{n, m, d} (B_{n,m} B_{m,n}) (ms_{n,d} ms_{d,n})
    if isinf(log_p)
        fprintf('Log_p is nan\n');
        keyboard;
    end
    if isnan(log_p)
        fprintf('Log_p is nan\n');
        keyboard;
    end
    log_p = -log_p;
    
    if nargout > 1
        
        % Prior Derivatives
        du = sum(mc, 2) - d_pts .* besseli(1, u, 1) ./ besseli(0, u, 1);
        % Chain rule to account for exponential
        dlog_u = du .* u;

        % Matrix derivatives
        dA = lambda2_y * (y * mc' - A * sum_ccT - B * sum_scT);
        
        dB = lambda2_y * (y * ms' - B * sum_ssT - A * sum_scT');

        % Precision Derivatives (used in chain rule)
        dlambda2_y = + 0.5 * d_pts * m_dim / lambda2_y ...
                     - 0.5 .* (trace(YYT) ...
                                - 2. * trace(mc' * A' * y) ...
                                - 2. * trace(ms' * B' * y) ...
                                + trace(ATA * sum_ccT) ...
                                + 2. * trace(ATB * sum_scT) ...
                                + trace(BTB * sum_ssT));

        % Chain rule to account for exponential
        dlog_lambda2_y = dlambda2_y * lambda2_y;
                
        dlog_p = -zm_pack_to_var_array(dlog_u, dA, dB, dlog_lambda2_y);
    end
end
