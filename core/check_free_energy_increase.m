function not_passed = check_free_energy_increase(fq, fq_old)
%CHECK_FREE_ENERGY_INCREASE Summary of this function goes here
%   Detailed explanation goes here
    if fq < fq_old && abs(fq - fq_old) > 1.E-2 * abs(fq_old)
        disp('!!! Free energy decreased !!!')
        disp(fq - fq_old)
        disp('Decreased by: ')
        disp(['Percent: ', num2str(abs(fq - fq_old) / abs(fq_old) * 100), '%'])
        disp(['Absolute: ', num2str(abs(fq - fq_old))])
        throw(MException('DO_VEM:FreeEnergyDecrease', ...
                         'Free energy decreased'));
    end
    not_passed = 0;
end

