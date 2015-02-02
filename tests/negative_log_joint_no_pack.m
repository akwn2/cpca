function [log_p, dlog_p] = negative_log_joint_no_pack(var_array, y, mc, ms, mc2, ms2, msc, D, M, N)
%NEGATIVE_LOG_JOINT Negative of the log-joint of for the M-step
%   Negative log joint and gradients for the probabilistic model of CPCA
    
    [u, A, B, log_lambda2_y] = vararray_unpack(var_array, D, M);
    
    % General pre-computable relations
    lambda2_y = exp(log_lambda2_y);
%     u = exp(log_u);
    u = u';
    U = repmat(u, [1, N]);
    YYT = y * y';
    
    ATA = A' * A;
    ATB = A' * B;
    BTB = B' * B;
    
    ccT = mc * mc' + diag(- diag(mc * mc') + sum(mc2, 2));
    scT = ms * mc' + diag(- diag(ms * mc') + sum(msc, 2));
    ssT = ms * ms' + diag(- diag(ms * ms') + sum(ms2, 2));
    
    % Gradient of the average log joint calculation
    log_p = N * (D - M / 2.) * log(2. * pi) ... %               2 pi constant from prior and likelihood
        - N * sum(u + log(besseli(0, u, 1))) ... %                      Prior normalizing constant
        + trace(U' * mc) ... %                                              Prior "energy"
        + 0.5 * N * M * log(lambda2_y) ... %                         Likelihood normalizing constant
        - 0.5 * lambda2_y * (trace(YYT) ... %                                   {m, d} y_{m,d} y_{m,d}
                                - 2. * trace(mc' * A' * y) ... %                {n, m, d} (mc_{d,n} A_{n,m}) y_{m,d}
                                - 2. * trace(ms' * B' * y) ... %                {n, m, d} (ms_{d,n} B_{n,m}) y_{m,d}
                                + trace(ATA * ccT) ... %                    {n, m, d} (A_{n,m} A_{m,n}) (mc_{n,d} mc_{d,n})
                                + 2 * trace(ATB * scT) ... %                {n, m, d} (A_{n,m} B_{m,n}) (ms_{n,d} mc_{d,n})
                                + trace(BTB * ssT)); %                      {n, m, d} (B_{n,m} B_{m,n}) (ms_{n,d} ms_{d,n})
    log_p = -log_p;
    
    if nargout > 1
        
        % Prior Derivatives
        du = sum(mc, 2) - N .* besseli(1, u, 1) ./ besseli(0, u, 1);
        % Chain rule to account for exponential
%         dlog_u = du .* u;

        % Matrix derivatives
        dA = lambda2_y * (y * mc' - A * ccT - B * scT);
        
        dB = lambda2_y * (y * ms' - B * ssT - A * scT');

        % Precision Derivatives (used in chain rule)
        dlambda2_y = + 0.5 * N * M / lambda2_y ...
                     - 0.5 .* (trace(YYT) ...
                                - 2. * trace(mc' * A' * y) ...
                                - 2. * trace(ms' * B' * y) ...
                                + trace(ATA * ccT) ...
                                + 2. * trace(ATB * scT) ...
                                + trace(BTB * ssT));

        % Chain rule to account for exponential
        dlog_lambda2_y = dlambda2_y * lambda2_y;
                
        dlog_p = - vararray_pack(du, dA, dB, dlog_lambda2_y);
    end
end
