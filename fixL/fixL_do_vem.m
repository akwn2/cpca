function [u, A, B, lambda2_y, run_stats] = fixL_do_vem(y, u, A, B, lambda2_y, k1, k2, m1, m2, y_held)
%DO_VEM Perform variational Expectation-Maximization on the CPCA model
%   Detailed explanation goes here
    converged = false;
    iter = 0;
    disp('Initializing run...')
    
    [m_dim, n_dim] = size(A);
    d_pts = size(y, 2);
    
    % Fix if no validation set given
    if nargin < 11
        % hold 20% of the data
        d_held = fix(.2 * d_pts);
        indexes = randperm(d_pts);
        y_held = y(:, indexes(1:d_held));
        
        % Fix all the indexing
        y = y(:, indexes(d_held:end));
        k1 = k1(:, indexes(d_held:end));
        k2 = k2(:, indexes(d_held:end));
        m1 = m1(:, indexes(d_held:end));
        m2 = m2(:, indexes(d_held:end));
        
        % Fix the dimension
        d_pts = size(y, 2);
    else
        d_held = size(y_held, 2);
    end
        
    old_model_pars = fixL_pack_model_pars(u, A, B);
    
    mf_pars = pack_mf_pars(k1, k2, m1, m2);
    
    [mc, ms, mc_sq, ms_sq, msc] = update_trig(mf_pars);
    
    pars = fixL_pack_pars(y, mc, ms, mc_sq, ms_sq, msc, n_dim, m_dim, d_pts, lambda2_y);
    
    fq = fixL_get_model_free_energy(old_model_pars, pars, mf_pars);
    
    score = fixL_get_validation_error(old_model_pars, pars, y_held, d_held);
    
    disp('Main algorithm')
    tic;
    hist_fq(1) = fq;
    hist_score(1) = score;
    hist_iter(1) = 0;
    hist_model_pars{1} = fixL_pack_to_var_array(u, A, B);
    while ~converged
        iter = iter + 1;
        disp(['Iteration: ', num2str(iter)]);
        
        % Main steps
        fprintf('\nE-step\n')
        [pars, mf_pars] = fixL_do_e_step(old_model_pars, pars, mf_pars, fq);
        
        fprintf('\nM-step\n')
        model_pars = fixL_do_m_step(old_model_pars, pars, 35);
        
        fprintf('\nFree energy\n')
        fq_old = fq;
        fq = fixL_get_model_free_energy(model_pars, pars, mf_pars);
        
        fprintf('\nValidation error\n')
        score_old = score;
        score = fixL_get_validation_error(model_pars, pars, y_held, d_held);
        
        % Iteration results
        celldisp(model_pars)
        
        % Overall convergence
%         check_free_energy_increase(fq, fq_old);
        converged = fixL_check_convergence(model_pars, old_model_pars, fq, fq_old, score, score_old);
        
        % Parameter iteration update
        old_model_pars = model_pars;
        hist_fq(iter + 1) = fq;
        hist_score(iter + 1) = score;
        hist_iter(iter + 1) = iter;
        
        [u, A, B] = fixL_unpack_model_pars(model_pars);
        hist_model_pars{iter + 1} = fixL_pack_to_var_array(u, A, B);
        
%         if mod(iter, 5) == 0
%             plot_within_iterations(model_pars, pars);
%             keyboard;
%         end
    end
    elapsed_time = toc;
    run_stats = {hist_fq, hist_iter, hist_score, hist_model_pars, elapsed_time};
    
%     % Output
%     figure
%     subplot(2,1,1)
%     plot(hist_iter, hist_fq);
%     title('Free energy by iteration');
%     subplot(2,1,2)
%     plot(hist_iter, hist_score);
%     title('Reconstruction score by iteration');
%     
%     plot_within_iterations(model_pars, pars);
    
    [u, A, B] = fixL_unpack_model_pars(model_pars);
end
