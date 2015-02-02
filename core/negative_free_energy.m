function [nfq, dnfq] = negative_free_energy(vararray, y, N, M, D)
%NEGATIVE_FREE_ENERGY Function that calculates the free energy for the
% CPCA model and its derivatives with respect to its parameters and the
% parameters of the mean field approximation.

    % Unpacking parameters:
    % Model parameters
    u = vararray(1:D);
    A = reshape(vararray(D + 1: D + D * M), [M, D]);
    B = reshape(vararray(D + D * M + 1: D + 2 * D * M), [M, D]);
    c = vararray(D + 2 * D * M + 1: D + 2 * D * M + 2);
    lambda2_y = vararray(D + 2 * D * M + 3);
    
    % Mean field terms
    k1 = reshape(vararray(D + 2 * D * M + 4:D + 2 * D * M + D * N + 3), [D, N]);
    k2 = reshape(vararray(D + 2 * D * M + D * N + 4:D + 2 * D * M + 2 * D * N + 3), [D, N]);
    m1 = reshape(vararray(D + 2 * D * M + 2 * D * N + 4 : D + 2 * D * M + 3 * D * N + 3), [D, N]);
    m2 = reshape(vararray(D + 2 * D * M + 3 * D * N + 4 : D + 2 * D * M + 4 * D * N + 3), [D, N]);
    
    % Get the expectation of trigonometric functions under the mean field
    [mc, ms, mc2, ms2, msc] = update_trig({k1, k2, m1, m2});
    
    % Negative Log Joint
    if nargout == 1
        nlogp = negative_log_joint(packet, y, mc, ms, mc2, ms2, msc, D, M, N); 
    else
        [nlogp, ndlogp] = negative_log_joint(packet, y, mc, ms, mc2, ms2, msc, D, M, N); 
    end
    
    % Negative Entropy
    nh = - get_gvm_entropy({k1, k2, m1, m2});
    
    % Negative Free Energy
    nfq = nlogp + nh;
    
    if nargout > 1
        
        % We proceed now to finding the derivatives of the free energy with
        % respect to the model parameters. These are found using the chain
        % rule in order to save computation on the moments, which are
        % assumed expensive to compute (e.g. found by a fine grid
        % integration).
        
        dlogp_mc = zeros(D, N);
        dlogp_ms = zeros(D, N);
        dlogp_mc2 = zeros(D, N);
        dlogp_ms2 = zeros(D, N);
        dlogp_msc = zeros(D, N);
        
        T0 = zeros(D, N);
        T1 = zeros(D, N);
        T2 = zeros(D, N);
        T3 = zeros(D, N);
        T4 = zeros(D, N);
        for nn = 1:N
            for dd = 1:D
        
                % 1. Find the derivatives of the log joint with respect to
                % the expected trigonometric functions mean field
        
                dlogp_mc(dd, :) = u(dd) - 0.5 * lambda2_y * A(:, dd)' * ...
                    (2 * (C - y) + A(:, 1:D ~= dd) * mc(1:D ~= dd, :) +...
                                   B(:, 1:D ~= dd) * ms(1:D ~= dd,:));

                dlogp_ms(dd, :) = - 0.5 * lambda2_y * B(:, dd)' * ...
                    (2 * (C - y) + B(:, 1:D ~= dd) * ms(1:D ~= dd, :) +...
                                   A(:, 1:D ~= dd) * mc(1:D ~= dd, :));

                dlogp_mc2(dd, :) = - 0.5 * lambda2_y * A(:, dd)' * A(:, dd);

                dlogp_ms2(dd, :) = - 0.5 * lambda2_y * B(:, dd)' * B(:, dd);

                dlogp_msc(dd, :) = - 0.5 * lambda2_y * A(:, dd)' * B(:, dd);

                
                % 2. Find the derivatives of the trigonometric moments that
                % will be used
        
                T0(dd, nn) = get_trigm(k1(dd, nn), k2(dd, nn), ...
                                       m1(dd, nn), m2(dd, nn));
                T1(dd, nn) = get_trigm(k1(dd, nn), k2(dd, nn), ...
                                       m1(dd, nn), m2(dd, nn));
                T2(dd, nn) = get_trigm(k1(dd, nn), k2(dd, nn), ...
                                       m1(dd, nn), m2(dd, nn));
                T3(dd, nn) = get_trigm(k1(dd, nn), k2(dd, nn), ...
                                       m1(dd, nn), m2(dd, nn));
                T4(dd, nn) = get_trigm(k1(dd, nn), k2(dd, nn), ...
                                       m1(dd, nn), m2(dd, nn));

                                   
                % 3. Find the derivative of the expected trigonometric
                % functions under the mean field with respect to the
                % mean field parameters
                
                
                
            end
        end
        
        % Aggregate negative log joint derivatives
        ndlogp = cat(1, ndlogp, dlogp_k1, dlogp_k2, dlogp_m1, dlogp_m2);
        
        % Aggregate negative entropy derivatives
        ndh = - cat(1, zeros(D + 2 * D * M + 3, 1), ...
                    dh_k1, dh_k2, dh_m1, dh_m2);
        
        % Negative Free Energy Gradient
        dnfq = ndlogp + ndh;
    end
end