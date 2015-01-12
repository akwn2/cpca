function [mc, ms, mc_sq, ms_sq, msc] = update_trig(mf_pars)
%UPDATE_TRIG Updates the trigonometric moments of the mean-field GvMs
%   Updates the trigonometric moments of the mean field GvMs
    if any(any(isnan(mf_pars{1}))) || any(any(isnan(mf_pars{2}))) || any(any(isinf(mf_pars{1}))) || any(any(isinf(mf_pars{2})))
        keyboard;
    end

%     [mom0, mom1, mom2] = get_gvm_grid_moments(mf_pars);    
    [mom0, mom1, mom2] = get_gvm_series_vectorized(mf_pars);

    % Obtaining values
    mom1_re = real(mom1);
    mom2_re = real(mom2);
    mom1_im = imag(mom1);
    mom2_im = imag(mom2);

    % Trigonometric moments
    mc = mom1_re ./ mom0;
    ms = mom1_im ./ mom0;
    mc_sq = 0.5 + 0.5 .* mom2_re ./ mom0;
    ms_sq = 0.5 - 0.5 .* mom2_re ./ mom0;
    msc = 0.5 .* mom2_im ./ mom0;
    
    if any(any(isnan(mc))) || any(any(isnan(ms))) || any(any(isinf(mc))) || any(any(isinf(ms)))
        keyboard;
    end
end

