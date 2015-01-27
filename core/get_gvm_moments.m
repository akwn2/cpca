function [mom0, mom1, mom2] = get_gvm_moments(mf)
%GET_GVM_SERIES_VECTORIZED Summary of this function goes here
%   Detailed explanation goes here
    [k1, k2, m1, m2] = unpack_mf(mf);
    
    [D, N] = size(k1);
    k1 = reshape(k1, D * N, 1);
    k2 = reshape(k2, D * N, 1);
    m1 = reshape(m1, D * N, 1);
    m2 = reshape(m2, D * N, 1);
    
    
    % Determine the calculation type based on the value of k1 + k2
    sum_k1_k2 = k1 + k2;
    max_series = 0.; % setting this to zero makes all calculations by grid
    calc_quad = find(sum_k1_k2 > max_series);
    calc_series = find(sum_k1_k2 <= max_series);
    
    mom0 = zeros(size(k1));
    mom1 = zeros(size(k1));
    mom2 = zeros(size(k1));
    
    % Calculations by Quadrature
    %-----------------------------
    for ii = 1:length(calc_quad)
        % Zeroth moment (Partition function)
        f0 = @(x) exp(k1(calc_quad(ii)) .* cos(x - m1(calc_quad(ii))) ...
                      + k2(calc_quad(ii)) .* cos(2 .*(x - m2(calc_quad(ii)))) ...
                      - k1(calc_quad(ii)) - k2(calc_quad(ii)));

        mom0(calc_quad(ii)) = quadgk(f0, -pi, pi);

        % First moment
        f1 = @(x) exp(+1.j .* x + k1(calc_quad(ii)) .* cos(x - m1(calc_quad(ii))) ...
                                + k2(calc_quad(ii)) .* cos(2 .*(x - m2(calc_quad(ii)))) ...
                                - k1(calc_quad(ii)) - k2(calc_quad(ii)));

        mom1(calc_quad(ii)) = quadgk(f1, -pi, pi);

        % Second moment
        f2 = @(x) exp(+2.j .* x + k1(calc_quad(ii)) .* cos(x - m1(calc_quad(ii))) ...
                                + k2(calc_quad(ii)) .* cos(2 .*(x - m2(calc_quad(ii)))) ...
                                - k1(calc_quad(ii)) - k2(calc_quad(ii)));

        mom2(calc_quad(ii)) = quadgk(f2, -pi, pi);
    end

    % Calculations by Series
    %-----------------------------
    if ~isempty(calc_series)
        max_i = 25;
        idx = 1:max_i;
        
        k1s = k1(calc_series);
        k2s = k2(calc_series);
        m1s = m1(calc_series);
        m2s = m2(calc_series);
        
        IDX = repmat(idx, [length(calc_series), 1]);
        K1 = repmat(k1s, [1 max_i]);
        K2 = repmat(k2s, [1 max_i]);
        M1 = repmat(m1s, [1 max_i]);
        M2 = repmat(m2s, [1 max_i]);


        % Moment 0
        mom0(calc_series) = besseli(0, k1s, 1) .* besseli(0, k2s, 1) ...
            + 2. * sum(besseli(2 .* IDX, K1, 1) .* ...
                           cos(2. .* IDX .* (M1 - M2)) .*...
                           besseli(IDX, K2, 1), 2);

        mom0(calc_series) = mom0(calc_series) * 2. * pi;

        % Moment 1
        mom1(calc_series) = besseli(1, k1s, 1) .* besseli(0, k2s, 1) ...
            + sum((exp(+2.0i .* IDX .* (M1 - M2)) .* ...
                       besseli(2 .* IDX + 1, K1, 1) + ...
                   exp(-2.0i .* IDX .* (M1 - M2)) .*...
                       besseli(2 .* IDX - 1, K1, 1)) .* ...
                           besseli(IDX, K2, 1), 2);

        mom1(calc_series) = mom1(calc_series) .* 2. .* pi .* exp(1.0i .* m1s);

        % Moment 2
        mom2(calc_series) = besseli(2, k1s, 1) .* besseli(0, k2s, 1) +...
            besseli(0, k1s, 1) .* besseli(1, k2s, 1) .* ...
                exp(-2.0i * (m1s - m2s)) +...
            sum(exp(+2.0i .* IDX .* (M1 - M2)) .* ...
                    besseli(2 .* IDX + 2, K1, 1) .* besseli(IDX, K2, 1) + ...
                exp(-2.0j .* (IDX + 1) .* (M1 - M2)) .*...
                    besseli(2 .* IDX, K1, 1) .* besseli(IDX + 1, K2, 1), 2);

        mom2(calc_series) = mom2(calc_series) .* 2. .* pi .* exp(2.0i .* m1s);
    end
    
    % Just to be safe, catch the cases when even the grid fails.
    fix = find(mom0 <= 0);
    if ~isempty(fix)
        F = length(fix);
        
        fprintf('-- Warning: using crude moment approximation. --\n');
        fprintf('--          Use: %d out of %d distributions.  --\n', ...
                F, N * D);
        % For these pathological cases, crudely approximate by a sampling
        % scheme that approximates a GvM sampler.
        
        Q = 1E4; % number of bins
        S = 1E4; % number of samples we draw
        samples = zeros(F, S); % samples container
        
        theta = linspace(-pi, pi, Q);
        for ff = 1:F
            gvm(ff, :) = exp(k1(ff) .* cos(theta - m2(ff)) + ...
                             k2(ff) * cos(2 .* (theta - m2(ff))) ...
                             - k1(ff) - k2(ff));
            samples(ff,:) = randsample(theta, S, true, gvm(ff, :));
        end
        
        mom0(fix) = sum(gvm, 2);
        mom1(fix) = sum(cos(samples) + 1.j * sin(samples), 2);
        mom2(fix) = sum(cos(2 * samples) + 1.j * sin(2 * samples), 2);
    end
    
    % Reshape all objects to the appropriate size
    mom0 = reshape(mom0, [D, N]);
    mom1 = reshape(mom1, [D, N]);
    mom2 = reshape(mom2, [D, N]);
end
