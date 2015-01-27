function not_passed = check_free_energy_increase(fq, fq_old)
%CHECK_FREE_ENERGY_INCREASE Summary of this function goes here
%   Detailed explanation goes here
    if fq < fq_old && abs(fq - fq_old) > 1.E-2 * abs(fq_old)
        fprintf('!!! FREE ENERGY DECREASED !!!\n')
        fprintf('Diffence in magnitude: %1.4e\n', abs(fq - fq_old))
        keyboard;
    end
    not_passed = 0;
end