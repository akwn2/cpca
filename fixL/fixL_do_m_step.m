function model_pars = fixL_do_m_step(model_pars, pars, max_iter)
%DO_M_STEP Performs the m_step of the CPCA model
%   Detailed explanation goes here
    
    % Unpacking
    [u, A, B] = fixL_unpack_model_pars(model_pars);
    [~, ~, ~, ~, ~, ~, n_dim, m_dim, ~, ~] = fixL_unpack_pars(pars);
    
    % Reshaping for optimization
    var_array = fixL_pack_to_var_array(log(u), A, B);
    
    % Optimizing
    var_array = minimize(var_array, 'fixL_negative_log_joint', -max_iter, pars{:});
%     var_array = minimize(var_array, 'regularized_neg_log_p', -10^(1 + max_iter), pars{:});
    
    % Reshaping
    [new_u, new_A, new_B] = fixL_unpack_from_var_array(var_array, n_dim, m_dim);
    
    % Packing & output
    model_pars = fixL_pack_model_pars(exp(new_u), new_A, new_B);
end

