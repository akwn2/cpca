function model = do_m_step(model, pars, max_iter)
%DO_M_STEP Performs the m_step of the CPCA model
%   Detailed explanation goes here
    
    % Unpacking
    [u, A, B, c, lambda2_y] = unpack_model(model);
    [~, ~, ~, ~, ~, ~, D, M, ~] = unpack_pars(pars);
    
    % Reshaping for optimization
    vararray = pack_array(u, A, B, c, log(lambda2_y));
    
    % Optimizing
    vararray = minimize(vararray, 'negative_log_joint', - max_iter, pars{:});
    
    % Reshaping
    [u, A, B, c, log_lambda2_y] = unpack_array(vararray, D, M);
    
    % Packing & output
    model = pack_model(u, A, B, c, exp(log_lambda2_y));
end

