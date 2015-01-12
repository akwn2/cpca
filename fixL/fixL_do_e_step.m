function [pars, mf_pars] = fixL_do_e_step(model_pars, pars, mf_pars, fq)
%DO_E_STEP E-step of the variational EM for CPCA
%   Calculates E-step using the mean field appromodel_parsimation with GvMs
    [u, A, B] = fixL_unpack_model_pars(model_pars);
    [y, mc, ms, ~, ~, ~, n_dim, m_dim, d_pts, lambda2_y] = fixL_unpack_pars(pars);
    
    % random permtations
    perm_hidden = randperm(n_dim);  
    [k1, k2, m1, m2] = unpack_mf_pars(mf_pars);

    for nn = 1:n_dim
        kk = perm_hidden(nn);

        rhs1 = u(kk) - lambda2_y * A(:, kk)' * ...
            (- y + A(:, 1:n_dim ~= kk) * mc(1:n_dim ~= kk, :) + B(:, 1:n_dim ~= kk) * ms(1:n_dim ~= kk,:));

        rhs2 = - lambda2_y * B(:, kk)' * ...
            (- y + B(:, 1:n_dim ~= kk) * ms(1:n_dim ~= kk, :) + A(:, 1:n_dim ~= kk) * mc(1:n_dim ~= kk,:));

        rhs3 = - 0.25 * lambda2_y * (A(:, kk)' * A(:, kk) - B(:, kk)' * B(:, kk));

        rhs4 = - 0.50 * lambda2_y * A(:, kk)' * B(:, kk);

        % Calculate the mean-field parameters

        z1 = rhs1 + 1.0j * rhs2;
        z2 = rhs3 + 1.0j * rhs4;

        k1(kk, :) = abs(z1);
        k2(kk, :) = abs(z2);
        m1(kk, :) = angle(z1);
        m2(kk, :) = 0.5 * angle(z2);

        mf_pars = pack_mf_pars(k1, k2, m1, m2);
        [mc, ms, mc_sq, ms_sq, msc] = update_trig(mf_pars);
        pars = fixL_pack_pars(y, mc, ms, mc_sq, ms_sq, msc, n_dim, m_dim, d_pts, lambda2_y);

        fq_old = fq;
        fq = fixL_get_model_free_energy(model_pars, pars, mf_pars);
%         check_free_energy_increase(fq, fq_old);
    end
end

