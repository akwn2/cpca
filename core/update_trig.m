function [mc, ms, mc2, ms2, msc] = update_trig(mf)
%UPDATE_TRIG Updates the trigonometric moments of the mean-field GvMs
%   Updates the trigonometric moments of the mean field GvMs
    if any(any(isnan(mf{1}))) || any(any(isnan(mf{2}))) || any(any(isinf(mf{1}))) || any(any(isinf(mf{2})))
        keyboard;
    end

%     [mom0, mom1, mom2] = get_gvm_grid_moments(mf);    
    [mom0, mom1, mom2] = get_gvm_moments(mf);

    % Obtaining values
    mom1_re = real(mom1);
    mom2_re = real(mom2);
    mom1_im = imag(mom1);
    mom2_im = imag(mom2);

    % Trigonometric moments
    mc = mom1_re ./ mom0;
    ms = mom1_im ./ mom0;
    mc2 = 0.5 + 0.5 .* mom2_re ./ mom0;
    ms2 = 0.5 - 0.5 .* mom2_re ./ mom0;
    msc = 0.5 .* mom2_im ./ mom0;
    
    if any(any(isnan(mc))) || any(any(isnan(ms))) || any(any(isinf(mc))) || any(any(isinf(ms)))
        keyboard;
    end
end

