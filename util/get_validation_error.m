function score = get_validation_error(model_pars, pars, y_held, d_held)
%GET_VALIDATION ERROR Get RMS error ini the reconstruction
%   Detailed explanation goes here

    % Unpacking model pars and replacing the data with held data for pars 
    [~, ~, A, B, ~] = unpack_model_pars(model_pars);
    [~, ~, ~, ~, ~, ~, n_dim, m_dim, ~] = unpack_pars(pars);
    
    pars = pack_pars(y_held, rand(n_dim, d_held), rand(n_dim, d_held), ...
        rand(n_dim, d_held), rand(n_dim, d_held), rand(n_dim, d_held), ...
        n_dim, m_dim, d_held);
    
    mf_pars = pack_mf_pars(rand(n_dim, d_held), rand(n_dim, d_held),...
        rand(n_dim, d_held), rand(n_dim, d_held));
    
    % Do an E-step to obtain mean field estimates
    pars = do_e_step(model_pars, pars, mf_pars, -1E20);
    [~, mc, ms, ~, ~, ~, ~, ~, ~] = unpack_pars(pars);
    
    % Predict
    y_test = A * mc + B * ms;
    
    % Score
    score = sqrt(sum(sum((y_held - y_test).^2)) / (m_dim * d_held));
    fprintf('Score: %3.3e\n', score);
end