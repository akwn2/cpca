function [log_p, dlog_p] = negative_log_joint_semi_vec(vararray, var_list, rows, cols, y, mc, ms, mc2, ms2,...
                                            msc, u, v, A, B, c, lambda2_y)
%NEGATIVE_LOG_JOINT_SEMI_VEC A semi-vectorised negative log-joint function for
% comparing with the results from the negative log joint to be used in the
% EM algorithm.

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
    
    [M, N] = size(y);
    D = size(A, 2);
    
    k = abs(u + 1.j * v);
    C = zeros(M, 1);
    C(1:2:end) = c(1);
    C(2:2:end) = c(2);

    prior = 0.;
    normTerm = 0.;
    
    for nn = 1:N
        % prior terms calculation
        prior = prior + u' * mc(:, nn) + v' * ms(:, nn);

        % likelihood terms calculation
        normTerm = normTerm + y(:, nn)' * y(:, nn) ... % YTY terms
            - 2 * y(:, nn)' * (A * mc(:, nn) + B * ms(:, nn) + C) ... % YT(model) terms
            + mc(:, nn)' * (A' * A - diag(diag(A' * A))) * mc(:, nn) ... % ATA off-diagonal terms
            + mc2(:, nn)' * diag(A' * A) ... % ATA diagonal terms
            + 2 * mc(:, nn)' * (A' * B - diag(diag(A' * B))) * ms(:, nn) ... % ATB off-diagonal terms
            + 2 * msc(:, nn)' * diag(A' * B) ... % ATB diagonal terms
            + ms(:, nn)' * (B' * B - diag(diag(B' * B))) * ms(:, nn) ... % BTB off-diagonal terms
            + ms2(:, nn)' * diag(B' * B) ... % BTB diagonal terms
            + (2 * A * mc(:, nn) + 2 * B * ms(:, nn) + C)' * C; % C terms
    end
    
    % log p calculation
    log_p = N * (- sum(k + log(2 * pi * besseli(0, k, 1))) ...
                   + M / 2 * log(lambda2_y / (2 * pi)))...
            + prior - 0.5 * lambda2_y * normTerm;

    log_p = - log_p;

    if nargout > 1
        du = zeros(size(u));
        dv = zeros(size(v));
        dA = zeros(size(A));
        dB = zeros(size(B));
        dC = zeros(size(C));
        dc = zeros(size(c));
        dlambda2_y = zeros(size(lambda2_y));
        
        for nn = 1:N
            du = du + mc(:, nn) ...
                - besseli(1, k, 1) ./ besseli(0, k, 1) .* u ./ k;
            
            dv = dv + ms(:, nn) ...
                - besseli(1, k, 1) ./ besseli(0, k, 1) .* v ./ k;
            
            dA = dA - lambda2_y .* ((C - y(:, nn)) * mc(:, nn)' ...
                                    + A * (mc(:, nn) * mc(:, nn)' ...
                                            - diag(diag(mc(:,nn) * mc(:, nn)'))...
                                            + diag(mc2(:, nn))) ...
                                    + B * (ms(:, nn) * mc(:, nn)' ...
                                            - diag(diag(ms(:, nn) * mc(:, nn)'))...
                                            + diag(msc(:, nn))));
            
            
            dB = dB - lambda2_y .* ((C - y(:, nn)) * ms(:, nn)' ...
                                    + B * (ms(:, nn) * ms(:, nn)' ...
                                            - diag(diag(ms(:,nn) * ms(:, nn)'))...
                                            + diag(ms2(:, nn))) ...
                                    + A * (mc(:, nn) * ms(:, nn)' ...
                                            - diag(diag(mc(:, nn) * ms(:, nn)'))...
                                            + diag(msc(:, nn))));
            
            dC = dC - lambda2_y .* (- y(:, nn) + (A * mc(:, nn) + B * ms(:, nn)) + C);
        end
        dc(1) = sum(dC(1:2:end));
        dc(2) = sum(dC(2:2:end));
        
        dlambda2_y = N * M / 2 * 1 / lambda2_y - 0.5 * normTerm;
        dlog_lambda2_y = dlambda2_y * lambda2_y;
        
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