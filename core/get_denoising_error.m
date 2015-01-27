function score = get_denoising_error(model, pars, y_held)
%GET_DENOISING_ERROR Get denoising error from the reconstruction after
%doing a E-step

    % Repackaging everything
    [~, A, B, c, ~] = unpack_model(model);
    [~, ~, ~, ~, ~, ~, D, M, ~] = unpack_pars(pars);
    
    N_held = size(y_held, 2);
    C = repmat(c, [M/2, N_held]);
    
    pars = pack_pars(y_held, rand(D, N_held), rand(D, N_held), ...
                     rand(D, N_held), rand(D, N_held),rand(D, N_held), ...
                     D, M, N_held);
    
    mf = pack_mf(rand(D, N_held), rand(D, N_held),...
                 rand(D, N_held), rand(D, N_held));
    
    % Obtain trigonometric moments from E-step for the prediction
    pars = do_e_step(model, pars, mf, -1E20, 2, false);
    [~, mc, ms, ~, ~, ~, ~, ~, ~] = unpack_pars(pars);
    
    % Predict the output
    y_pred = A * mc + B * ms + C;
    
    % Calculate RMS
    score = sqrt(norm(y_held - y_pred) ^ 2 / (M * N_held));
    fprintf('\t\t\tDenoising RMS: %3.3e\n', score);
end