function [pars, mf, fq] = do_e_step(model_pars, pars, ...
                                         mf, fq, eStepPasses, checkFq)
%DO_E_STEP E-step of the variational EM for CPCA
%   Calculates E-step using the mean field approximation with GvMs
    [u, A, B, c, lambda2_y] = unpack_model(model_pars);
    [y, mc, ms, ~, ~, ~, D, M, N] = unpack_pars(pars);
    
    C = repmat(c, [M / 2, N]);
    [k1, k2, m1, m2] = unpack_mf(mf);
    
    % We permute the order of the updates to avoid problems in this partial
    % e-step.
    perm_hidden = randperm(D);
    
    for ii = 1:eStepPasses
        fprintf('\t\tE-step pass: %d of %d.\n', ii, eStepPasses);

        for dd = 1:D
            kk = perm_hidden(dd);

            % Updates
            rhs1 = u(kk) - 0.5 * lambda2_y * A(:, kk)' * ...
                (2 * (C - y) + 2 * A(:, 1:D ~= kk) * mc(1:D ~= kk, :) +...
                               2 * B(:, 1:D ~= kk) * ms(1:D ~= kk,:));

            rhs2 = - 0.5 * lambda2_y * B(:, kk)' * ...
                (2 * (C - y) + 2 * B(:, 1:D ~= kk) * ms(1:D ~= kk, :) +...
                               2 * A(:, 1:D ~= kk) * mc(1:D ~= kk,:));

            rhs3 = - 0.25 * lambda2_y * (A(:, kk)' * A(:, kk) - ...
                                        B(:, kk)' * B(:, kk));

            rhs4 = - 0.5 * lambda2_y * A(:, kk)' * B(:, kk);

            % Calculate the mean-field parameters
            z1 = rhs1 + 1.0j * rhs2;
            z2 = rhs3 + 1.0j * rhs4;
            k1(kk, :) = abs(z1);
            k2(kk, :) = abs(z2);
            m1(kk, :) = angle(z1);
            m2(kk, :) = 0.5 * angle(z2);

            % Update moments
            mf = pack_mf(k1, k2, m1, m2);
            [mc, ms, mc_sq, ms_sq, msc] = update_trig(mf);
            pars = pack_pars(y, mc, ms, mc_sq, ms_sq, msc, D, M, N);
        end
    end
    % Update free energy
    if checkFq
        fq = get_model_free_energy(model_pars, pars, mf, fq);
    end
end

