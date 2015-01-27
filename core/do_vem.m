function [u, A, B, c, lambda2_y, run_stats] = do_vem(y, u, A, B, c, lambda2_y, k1, k2, m1, m2, y_held)
%DO_VEM Perform variational Expectation-Maximization on the CPCA model
%   Detailed explanation goes here
    
    % Main options
    fprintf('-----------------------------------------------\n');
    fprintf('Running variational EM algorithm for CPCA model\n');
    fprintf('-----------------------------------------------\n');
    eStepPasses = 5;
    rtol = 1.E-4;
    fq = -1E100;
    
    % Dimension of a single data point, hidden space and number of points
    [M, D] = size(A);
    N = size(y, 2);
    
    % Priming variables
    fprintf('Priming variables\n');
    model = pack_model(u, A, B, c, lambda2_y);
    mf = pack_mf(k1, k2, m1, m2);
    [mc, ms, mc2, ms2, msc] = update_trig(mf);
    pars = pack_pars(y, mc, ms, mc2, ms2, msc, D, M, N);
    fq = get_model_free_energy(model, pars, mf, fq);
    score = get_denoising_error(model, pars, y_held);
    
    
    % Main EM steps
    fprintf('Main algorithm\n');
    tic;
    iter = 1;
    converged = false;
    while ~converged
        fprintf('\tIteration: %d\n', iter);
        
        fprintf('\tE-step\n')
        [pars, mf, fq] = do_e_step(model, pars,...
                                   mf, fq, eStepPasses, true);
        
        fprintf('\tM-step\n')
        old_model = model;
        model = do_m_step(old_model, pars, 30);
        
        fprintf('\tFree energy\n')
        [fq, fq_old] = get_model_free_energy(model, pars, mf, fq);
        
        fprintf('\tPredict (do a denoising)\n')
        score_old = score;
        score = get_denoising_error(model, pars, y_held);
        
        % Display iteration results
        celldisp(model)
        
        % Check convergence in the model parameters, free energy and
        % prediction error
        converged = check_convergence(rtol, model, old_model,...
                                      fq, fq_old, score, score_old);
        
        % Keep history of key quantities for later analysis
        hist_fq(iter) = fq;
        hist_denoise(iter) = score;
        hist_model{iter} = model;
        iter = iter + 1;
    end
    
    elapsed_time = toc;
    run_stats = {hist_fq, hist_denoise, hist_model, elapsed_time};
    [u, A, B, c, lambda2_y] = unpack_model(model);
end