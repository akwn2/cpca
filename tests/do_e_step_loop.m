function [mc, ms, mc2, msc, ms2, k1, k2, m1, m2, fq, fq_old] = do_e_step_loop(y, ...
                                              u, v, A, B, c, lambda2_y, ...
                                              mc, ms, mc2, ms2, msc, ...
                                              k1, k2, m1, m2, fq, fq_old,...
                                              eStepPasses, var_list, checkFq)
%DO_E_STEP_LOOP E-step of the variational EM for CPCA using loops
%   Calculates E-step using the mean field appromodel_parsimation with GvMs

    [M, N] = size(y);
    D = size(A,2);
    
    C = zeros(M,1);
    C(1:2:end) = c(1);
    C(2:2:end) = c(2);
    
    % random permutations
%     perm_hidden = randperm(D);
     for ii = 1:eStepPasses
        fprintf('\t\tE-step pass: %d of %d.\n', ii, eStepPasses);
        for nn = 1:N
            for dd = 1:D
                z1c = 0;
                z1s = 0;
                z2c = 0;
                z2s = 0;
                
                for mm = 1:M
                    w_sum = 0.;
                    for ww = 1:D
                        if ww ~= dd
                            w_sum = w_sum + A(mm, ww) * mc(ww, nn) + ...
                                            B(mm, ww) * ms(ww, nn);
                        end
                    end
                    
                    z1c = z1c + A(mm, dd) * (2 * (C(mm) - y(mm, nn)) + w_sum);
                    z1s = z1s + B(mm, dd) * (2 * (C(mm) - y(mm, nn)) + w_sum);
                    z2c = z2c + A(mm, dd) ^ 2 - B(mm, dd) ^ 2;
                    z2s = z2s + A(mm, dd) * B(mm, dd); 
                end
                
                z1c = u(dd) - 0.5 * lambda2_y * z1c;
                z1s = v(dd) - 0.5 * lambda2_y * z1s;
                z2c = - 0.5 * lambda2_y * z2c;
                z2s = - 0.5 * lambda2_y * z2s;
                
                z1 = z1c + 1.0j * z1s;
                z2 = z2c + 1.0j * z2s;                

                k1(dd, nn) = abs(z1);
                k2(dd, nn) = abs(z2);
                m1(dd, nn) = angle(z1);
                m2(dd, nn) = 0.5 * angle(z2);
                
                [mc, ms, mc2, ms2, msc] = update_trig(k1, k2, m1, m2);
            end
        end
        if checkFq
            [fq, fq_old] = get_model_free_energy(y, u, v, A, B, c, lambda2_y, ...
                                                mc, ms, mc2, ms2, msc, ...
                                                k1, k2, m1, m2, fq, ...
                                                fq_old, var_list);
        end
     end
end

