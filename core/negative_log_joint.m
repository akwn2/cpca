function [log_p, dlog_p] = negative_log_joint(vararray, y, mc, ms, mc2, ms2, msc, D, M, N)
%NEGATIVE_LOG_JOINT Negative of the log-joint of for the M-step
%   Negative log joint and gradients for the probabilistic model of CPCA
    
    [u, A, B, c, log_lambda2_y] = unpack_array(vararray, D, M);
    
    % General pre-computable relations
    C = repmat(c, [floor(M / 2), N]);
    lambda2_y = exp(log_lambda2_y);
    
    ATA = A' * A;
    ATB = A' * B;
    BTB = B' * B;
    
    % Prior terms
    prior = sum(u' * mc);
    
    % Norm within the Liuelihood term
    ccT = (mc * mc') + diag(- diag(mc * mc') + sum(mc2, 2));
    scT = (ms * mc') + diag(- diag(ms * mc') + sum(msc, 2));
    ssT = (ms * ms') + diag(- diag(ms * ms') + sum(ms2, 2));
    
    normTerm = trace((y - 2 * (A * mc + B * ms + C)) * y') ... % y terms
            + trace(ATA * ccT + 2 * ATB * scT + BTB * ssT) + ... % <cos'A'Acos> + <cos'A'Bsin> + <sin'B'Bsin>
            + trace((2 * (A * mc + B * ms) + C) * C');%        (A <cos> + B <sin> + C)' * C
        
    % Log joint calculation
    log_p = N * (-sum(u + log(2 * pi * besseli(0, u, 1))) ...
                   + M / 2 * log(lambda2_y / (2 * pi)))...
            + prior - 0.5 * lambda2_y * normTerm;
    
    log_p = -log_p;
    
    if nargout > 1
        
        % Prior Derivatives
        du = sum(mc, 2) - N .* besseli(1, u, 1) ./ besseli(0, u, 1);

        % Matrix derivatives
        dA = -lambda2_y * ((C - y) * mc' + A * ccT + B * scT);
        dB = -lambda2_y * ((C - y) * ms' + B * ssT + A * scT');
        
        % Offset derivative
        aux = -lambda2_y * sum(- y + A * mc + B * ms + C, 2);
        
        % shared c's for two-dimensional data
        dc = zeros(size(c));
        dc(1) = sum(aux(1:2:end));
        dc(2) = sum(aux(2:2:end));

        % Precision Derivatives (used in chain rule)
        dlambda2_y = + 0.5 * (N * M / lambda2_y - normTerm);

        % Chain rule to account for exponential
        dlog_lambda2_y = dlambda2_y * lambda2_y;
        
        dlog_p = -cat(1, du, reshape(dA, [M * D, 1]), ...
                    reshape(dB, [M * D, 1]), dc, dlog_lambda2_y);
    end
end