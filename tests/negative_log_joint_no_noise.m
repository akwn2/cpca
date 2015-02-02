function [log_p, dlog_p] = negative_log_joint_no_noise(vararray, var_list, rows, cols, y, mc, ms, mc2, ms2, msc, u, v, A, B, c, lambda2_y)
%NEGATIVE_LOG_JOINT_NO_NOISE Negative of the log-joint of for the M-step
%without the noise (fixed noise)
%   Negative log joint and gradients for the probabilistic model of CPCA
    [M, N] = size(y);
    D = size(A,2);
    
    % with this unpacking we overwrite all variables and leave the ones we
    % want to keep fixed as the value supplied as a parameter.
    unpacked = unpackit(vararray, rows, cols);
    for ii = 1:length(var_list)
        % account for the case you have positivity constraint
        if strcmp(var_list{ii}, 'lambda2_y')
            eval([var_list{ii}, '= exp(unpacked{ii});']);
        else
            eval([var_list{ii}, '= unpacked{ii};']);
        end
    end
    
    % General pre-computable relations
%     u = exp(log_u);
    
    U = repmat(u, [1, N]);
    YYT = y * y';
    
    ATA = A' * A;
    ATB = A' * B;
    BTB = B' * B;
    
    sum_ccT = mc * mc' + diag(- diag(mc * mc') + sum(mc2, 2));
    sum_scT = ms * mc' + diag(- diag(ms * mc') + sum(msc, 2));
    sum_ssT = ms * ms' + diag(- diag(ms * ms') + sum(ms2, 2));
    
    % Gradient of the average log joint calculation
    log_p = N * (D - M / 2.) * log(2. * pi) ... %               2 pi constant from prior and likelihood
        - N * sum(u + log(besseli(0, u, 1))) ... %                      Prior normalizing constant
        + trace(U' * mc) ... %                                              Prior "energy"
        + 0.5 * N * M * log(lambda2_y) ... %                         Likelihood normalizing constant
        - 0.5 * lambda2_y * (trace(YYT) ... %                                   sum_{m, d} y_{m,d} y_{m,d}
                                - 2. * trace(mc' * A' * y) ... %                sum_{n, m, d} (mc_{d,n} A_{n,m}) y_{m,d}
                                - 2. * trace(ms' * B' * y) ... %                sum_{n, m, d} (ms_{d,n} B_{n,m}) y_{m,d}
                                + trace(ATA * sum_ccT) ... %                    sum_{n, m, d} (A_{n,m} A_{m,n}) (mc_{n,d} mc_{d,n})
                                + 2 * trace(ATB * sum_scT) ... %                sum_{n, m, d} (A_{n,m} B_{m,n}) (ms_{n,d} mc_{d,n})
                                + trace(BTB * sum_ssT)); %                      sum_{n, m, d} (B_{n,m} B_{m,n}) (ms_{n,d} ms_{d,n})

    log_p = -log_p;
    
    if nargout > 1
        
        % Prior Derivatives
        du = sum(mc, 2) - N .* besseli(1, u, 1) ./ besseli(0, u, 1);
        % Chain rule to account for exponential
%         dlog_u = du .* u;

        % Matrix derivatives
        dA = lambda2_y * (y * mc' - A * sum_ccT - B * sum_scT);
        
        dB = lambda2_y * (y * ms' - B * sum_ssT - A * sum_scT');
                

        % Create a string to evaluate the correct gradients
        eval_str = ['{'];
        for ii = 1:length(var_list)
            % account for the case you have positivity constraint
            if strcmp(var_list{ii}, 'lambda2_y')
                eval_str = [eval_str, ', dlog_', var_list{ii}];
            else
                eval_str = [eval_str, ' , d', var_list{ii}];
            end
        end
        eval_str = [eval_str, '}'];
        to_pack = eval(eval_str);
        
        dlog_p = - packit(to_pack{:})';
    end
end
