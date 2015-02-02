function [mc, ms, mc2, msc, ms2, k1, k2, m1, m2, fq, fq_old] = do_e_step_semi_vec(y, ...
                                              u, v, A, B, c, lambda2_y, ...
                                              mc, ms, mc2, ms2, msc, ...
                                              k1, k2, m1, m2, fq, fq_old,...
                                              eStepPasses, var_list, checkFq)
%DO_E_STEP_LOOP E-step of the variational EM for CPCA using loops
%   Calculates E-step using the mean field appromodel_parsimation with GvMs

    [M, N] = size(y);
    D = size(A,2);
    
    C = zeros(M,N);
    C(1:2:end,:) = c(1);
    C(2:2:end,:) = c(2);
    
    % random permutations
%     perm_hidden = randperm(D);
     for ii = 1:eStepPasses
        fprintf('\t\tE-step pass: %d of %d.\n', ii, eStepPasses);
        for dd = 1:D
            w_sum = zeros(M, N);
            for ww = 1:D
                if ww ~= dd
                    w_sum = w_sum + A(:, ww) * mc(ww, :) + ...
                                    B(:, ww) * ms(ww, :);
                end
            end

            z1c = A(:, dd)' * (2 * (C - y) + w_sum);
            z1s = B(:, dd)' * (2 * (C - y) + w_sum);
            z2c = sum(A(:, dd) .^ 2 - B(:, dd) .^ 2);
            z2s = sum(A(:, dd) .* B(:, dd)); 

            z1c = u(dd) * ones(1, N) - 0.5 * lambda2_y * z1c;
            z1s = v(dd) * ones(1, N) - 0.5 * lambda2_y * z1s;
            z2c = - 0.5 * lambda2_y * z2c * ones(1, N);
            z2s = - 0.5 * lambda2_y * z2s * ones(1, N);

            z1 = z1c + 1.0j * z1s;
            z2 = z2c + 1.0j * z2s;                

            k1(dd, :) = abs(z1);
            k2(dd, :) = abs(z2);
            m1(dd, :) = angle(z1);
            m2(dd, :) = 0.5 * angle(z2);

            [mc, ms, mc2, ms2, msc] = update_trig(k1, k2, m1, m2);
        end
        if checkFq
            [fq, fq_old] = get_model_free_energy(y, u, v, A, B, c, lambda2_y, ...
                                                mc, ms, mc2, ms2, msc, ...
                                                k1, k2, m1, m2, fq, ...
                                                fq_old, var_list);
        end
     end
end

