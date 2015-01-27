function h = get_gvm_entropy(mf)
%GET_GVM_ENTROPY Calculates the differential entropy for the MF GvM's
%   Detailed explanation goes here
    [k1, k2, m1, m2] = unpack_mf(mf);
    
    % Get moments (exponentially weighted)
    [mom0, mom1, mom2] = get_gvm_moments(mf);
    mom1_re = real(mom1);
    mom2_re = real(mom2);
    mom1_im = imag(mom1);
    mom2_im = imag(mom2);
    
    h = sum(sum(log(mom0) + k1 + k2 +...
            - (k1 .* (cos(m1) .* mom1_re  + sin(m1) .* mom1_im)...
               + k2 .* (cos(2. * m2) .* mom2_re + ...
                        sin(2. * m2) .* mom2_im)) ./ mom0));
    if isnan(h)
        keyboard;
    end
end

