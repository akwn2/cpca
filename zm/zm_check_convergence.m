function converged = zm_check_convergence(model_pars, old_model_pars, fq, fq_old, score, score_old)
%CHECK_CONVERGENCE Checks convergence of vem iterations.
%   Detailed explanation goes here
    [u, A, B, lambda_y] = zm_unpack_model_pars(model_pars);
    var_array = zm_pack_to_var_array(u, A, B, lambda_y);
    
    [u_old, A_old, B_old, lambda_y_old] = zm_unpack_model_pars(old_model_pars);
    old_var_array = zm_pack_to_var_array(u_old, A_old, B_old, lambda_y_old);
    
    model_res = norm(var_array - old_var_array);
    fq_res = norm(fq - fq_old);
    score_res = norm(score - score_old);
    
    converged = (model_res < 1.E-2 * norm(old_var_array)) || ...
                   (fq_res < 1.E-3 * abs(fq_old)) || ...
                   (score_res < 5.E-3 * abs(score_old));
end

