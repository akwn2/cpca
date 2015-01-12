function [mom0, mom1, mom2] = gvm_series_moments(mf_pars)
%GVM_SERIES_MOMENTS Calculates the GvM moments using series expansion formulae
%   Detailed explanation goes here
    [k1, k2, m1, m2] = unpack_mf_pars(mf_pars);

    [dim_n, d_pts] = size(k1);
    mom0 = zeros(dim_n, d_pts);
    mom1 = zeros(dim_n, d_pts);
    mom2 = zeros(dim_n, d_pts);
    
    for nn = 1:dim_n
        mom0(nn, :) = besseli(0, k1(nn, :), 1) .* besseli(0, k2(nn, :), 1);
        mom1(nn, :) = besseli(1, k1(nn, :), 1) .* besseli(0, k2(nn, :), 1);
        mom2(nn, :) = besseli(2, k1(nn, :), 1) .* besseli(0, k2(nn, :), 1) ...
            + besseli(0, k1(nn, :), 1) .* besseli(1, k2(nn, :), 1) .* exp(-2.0i * (m1(nn, :) - m2(nn, :)));
        
        max_i = min(10 + ceil(0.5 * (k1(nn, :) + k2(nn, :))), 50);
        
        [idx, kappa1] = meshgrid(1:max_i, k1(nn,:), 1:max_i);
        [~, kappa2] = meshgrid(1:max_i, k2(nn,:), 1:max_i);
        [~, delta] = meshgrid(1:max_i, m1(nn,:) - m2(nn,:));

        mom0(nn, :) = mom0(nn, :) + ...
            2. * sum(besseli(2.0 .* idx, kappa1) .* cos(2.0 .* idx .* delta) .* besseli(idx, kappa2, 1), 2)';

        mom1(nn, :) = mom1(nn, :) + sum((exp(+2.0i .* idx .* delta) .* besseli(2 .* idx + 1, kappa1) +...
                                         exp(-2.0i .* idx .* delta) .* besseli(2 .* idx - 1, kappa1)) .* besseli(idx, kappa2, 1), 2)';

        mom2(nn, :) = mom2(nn, :) + sum(exp(+2.0i .* idx .* delta) .* besseli(2.0 .* idx + 2, kappa1, 1) .* besseli(idx, kappa2, 1) +...
                                        exp(-2.0i .* (idx + 1) .* delta) .* besseli(2.0 .* idx, kappa1, 1) .* besseli(idx + 1, kappa2, 1), 2)';
                    
        % Moment without the exponential re-weighting
        mom0(nn, :) = 2. * pi * mom0(nn, :);
        mom1(nn, :) = 2. * pi * mom1(nn, :) .* exp(1.0i .* m1(nn, :));
        mom2(nn, :) = 2. * pi * mom2(nn, :) .* exp(2.0i .* m1(nn, :));
    end
end

