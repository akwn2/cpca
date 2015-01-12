function model_pars = init_sequential_do_m_step(model_pars, pars, max_iter)
%DO_M_STEP Performs the m_step of the CPCA model
%   Detailed explanation goes here
    
    % Unpacking
    [u, v, A, B, lambda2_y] = unpack_model_pars(model_pars);
    [~, ~, ~, ~, ~, ~, n_dim, m_dim, ~] = unpack_pars(pars);
    
    % Reshaping for optimization
    ufixed = u(1:end - 1);
    vfixed = v(1:end - 1);
    Afixed = A(:,1:end - 1);
    Bfixed = B(:,1:end - 1);
    
    uvar = u(end);
    vvar = v(end);
    Avar = A(:, end);
    Bvar = B(:, end);
    
    pars{10} = ufixed;
    pars{11} = vfixed;
    pars{12} = Afixed;
    pars{13} = Bfixed;
    
    var_array = init_sequential_pack_to_var_array(uvar, vvar, Avar, Bvar, log(lambda2_y));
    
    % Optimizing
    var_array = minimize(var_array, 'init_sequential_negative_log_joint', - max_iter, pars{:});
    
    % Reshaping
    [new_u, new_v, new_A, new_B, new_log_lambda2_y] = ...
        init_sequential_unpack_from_var_array(var_array, m_dim);
    
    % Packing & output
    new_u = [ufixed, new_u]';
    new_v = [vfixed, new_v]';
    new_A = [Afixed, new_A];
    new_B = [Bfixed, new_B];
    
    model_pars = pack_model_pars(new_u, new_v, new_A, new_B, exp(new_log_lambda2_y));
end

