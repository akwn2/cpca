function converged = check_convergence(rtol, model, old_model, fq, fq_old, score, score_old)
%CHECK_CONVERGENCE Checks convergence of vem iterations.
%   Detailed explanation goes here
    
    % Obtaining arrays to compute the norm of the difference in the model
    [u, v, A, B, lambda_y] = unpack_model(model);
    vararray = pack_array(u, v, A, B, lambda_y);
    
    [u_old, v_old, A_old, B_old, lambda_y_old] = unpack_model(old_model);
    old_vararray = pack_array(u_old, v_old, A_old, B_old, lambda_y_old);
    
    % Computing the convergence criterions
    model_res = norm(vararray - old_vararray);
    fq_res = norm(fq - fq_old);
    score_res = norm(score - score_old);
    
    converged = false;
    if model_res < rtol * norm(old_vararray)
        converged = true;
        fprintf('Converged in parameters\n');
    end
    if fq_res < rtol * abs(fq_old)
        converged = true;
        fprintf('Converged in free energy\n');
    end
    if score_res < rtol * abs(score_old)
        converged = true;
        fprintf('Converged in denoising\n');
    end
end

