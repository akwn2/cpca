function [fq, h] = get_model_free_energy(model, pars, mf, fq)
%GET_MODEL_FREE_ENERGY Calculates the model Free Energy
%   Detailed explanation goes here
    [u, A, B, c, lambda2_y] = unpack_model(model);
    var_array = pack_array(u, A, B, c, log(lambda2_y));
    
    fprintf('\t\tUpdating the free energy.\n');
    
    fq_old = fq;
    
    log_p = - negative_log_joint(var_array, pars{:});
    fprintf('\t\t\tlog p = %3.3e\n', log_p);
    
    h = get_gvm_entropy(mf);    
    fprintf('\t\t\th = %3.3e\n', h);
    
    fq = log_p + h;
    fprintf('\t\t\tfq = %3.3e\n', fq);
    
    if fq < fq_old
        fprintf('!!! FREE ENERGY DECREASED !!!\n')
        fprintf('Old Fq - new Fq = %1.4e\n', fq_old - fq)
        keyboard;
    end
end

