function [log_p, dlog_p] = negative_log_joint_fd(vararray, var_list, rows, cols, y, mc, ms, mc2, ms2, msc, u, v, A, B, c, lambda2_y)
%NEGATIVE_LOG_JOINT_FD Negative of the log-joint of for the M-step using
%finite differences to find the gradients.
%   Negative log joint and gradients for the probabilistic model of CPCA
    
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
    
    [M, ~] = size(A);
    N = size(y,2);
    
    % General pre-computable relations
    C = repmat(c, [floor(M / 2), N]);
    k = abs(u + 1.i * v);
    
    ATA = A' * A;
    ATB = A' * B;
    BTB = B' * B;
    
    sum_mc2 = sum(mc2, 2);
    sum_msc = sum(msc, 2);
    sum_ms2 = sum(ms2, 2);
    
    % Prior terms
    prior = sum(u' * mc + v' * ms);
    
    % Norm within the Likelihood term
    ccT = (mc * mc') + diag(- diag(mc * mc') + sum_mc2);
    scT = (ms * mc') + diag(- diag(ms * mc') + sum_msc);
    ssT = (ms * ms') + diag(- diag(ms * ms') + sum_ms2);
    
    normTerm = trace((y - 2 * (A * mc + B * ms + C)) * y') ... % y terms
            + trace(ATA * ccT + 2 * ATB * scT + BTB * ssT) + ... % <cos'A'Acos> + <cos'A'Bsin> + <sin'B'Bsin>
            + trace((2 * (A * mc + B * ms) + C) * C');%        (A <cos> + B <sin> + C)' * C
        
    % Log joint calculation
    log_p = N * (-sum(k + log(2 * pi * besseli(0, k, 1))) ...
                   + M / 2 * log(lambda2_y / (2 * pi)))...
            + prior - 0.5 * lambda2_y * normTerm;
    
    log_p = -log_p;
    
    if nargout > 1
        
        % Find the derivatives using finite differences
        e = 1.E-8;
        dlog_p = zeros(length(vararray),1);
        
        for j = 1:length(vararray)
          dx = zeros(size(vararray));
          dx(j) = dx(j) + e; % perturbation
          % Forward point
          y2 = negative_log_joint_fd(vararray + dx,...
                                      var_list, rows, cols, y,...
                                      mc, ms, mc2, ms2, msc,...
                                      u, v, A, B, c, lambda2_y);
          dx = -dx;
          % Backward point
          y1 = negative_log_joint_fd(vararray + dx,...
                                      var_list, rows, cols, y,...
                                      mc, ms, mc2, ms2, msc,...
                                      u, v, A, B, c, lambda2_y);
          dlog_p(j) = (y2 - y1)/(2*e);
        end
        
    end
end