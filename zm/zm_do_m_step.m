function model_pars = zm_do_m_step(model_pars, pars, max_iter)
%DO_M_STEP Performs the m_step of the CPCA model
%   Detailed explanation goes here
    
    % Unpacking
    [u, A, B, lambda2_y] = zm_unpack_model_pars(model_pars);
    [~, ~, ~, ~, ~, ~, n_dim, m_dim, ~] = unpack_pars(pars);
    
    % Reshaping for optimization
    var_array = zm_pack_to_var_array(log(u), A, B, log(lambda2_y));
    
    % Optimizing
    var_array = minimize(var_array, 'zm_negative_log_joint', -max_iter, pars{:});
%     var_array = minimize(var_array, 'regularized_neg_log_p', -10^(1 + max_iter), pars{:});
    
    % Reshaping
    [new_u, new_A, new_B, new_log_lambda2_y] = zm_unpack_from_var_array(var_array, n_dim, m_dim);
    
    if new_log_lambda2_y < log(0.025)
        new_log_lambda2_y = log(0.025);
    end
    % Packing & output
    model_pars = zm_pack_model_pars(exp(new_u), new_A, new_B, exp(new_log_lambda2_y));
end

