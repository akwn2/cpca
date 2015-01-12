function fq = zm_get_model_free_energy(model_pars, pars, mf_pars)
%GET_MODEL_FREE_ENERGY Calculates the model Gibbs Energy
%   Detailed explanation goes here
    [u, A, B, lambda2_y] = zm_unpack_model_pars(model_pars);
    var_array = zm_pack_to_var_array(log(u), A, B, log(lambda2_y));
    
    log_p = - zm_negative_log_joint(var_array, pars{:});
    fprintf('log p = %3.3e\n', log_p);
    
    h = get_gvm_entropy(mf_pars);
    fprintf('h = %3.3e\n', h);
    
    fq = log_p + h;
    fprintf('fq = %3.3e\n', fq);
end

