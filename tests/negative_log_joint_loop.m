function log_p = negative_log_joint_loop(y, mc, ms, mc2, ms2,...
                                            msc, u, v, A, B, c, lambda2_y)
%NEGATIVE_LOG_JOINT_LOOP A non-vectorised negative log-joint function for
% comparing with the results from the negative log joint to be used in the
% EM algorithm.
    [M, N] = size(y);
    D = size(A, 2);
    
    k = abs(u + 1.j * v);
    C = zeros(M, 1);
    C(1:2:end) = c(1);
    C(2:2:end) = c(2);

    % prior terms calculation
    prior = 0.;
    for nn = 1: N
        for dd = 1:D
            prior = prior + u(dd) * mc(dd, nn) + v(dd) * ms(dd, nn);
        end
    end

    % likelihood terms calculation
    normTerm = 0;
    for nn = 1: N
        for mm = 1:M
            sum_over_pp = 0;
            sum_over_ij = 0;
            sum_over_qq = 0;
            sum_over_zz = 0;
            
            for dd = 1:D
                % Sum of the model predictions without C, i.e., A*mc+B*ms
                sum_over_pp = sum_over_pp + A(mm, dd) * mc(dd, nn) + ...
                                            B(mm, dd) * ms(dd, nn);

                % Sum of the off-diagonal terms on the quadratic
                for jj = 1:D
                    if jj ~= dd
                        sum_over_ij = sum_over_ij + ...
                               + mc(dd, nn) * A(mm, dd) * A(mm, jj) * mc(jj, nn) ...
                            +2 * mc(dd, nn) * A(mm, dd) * B(mm, jj) * ms(jj, nn) ...
                               + ms(dd, nn) * B(mm, dd) * B(mm, jj) * ms(jj, nn);
                    end 
                end
                
                % Sum of the diagonal terms in the quadratics
                sum_over_qq = sum_over_qq + A(mm, dd) ^ 2 * mc2(dd, nn) ...
                                          + B(mm, dd) ^ 2 * ms2(dd, nn) ...
                                          + 2 * A(mm, dd) * B(mm, dd) * msc(dd, nn);

                % Sum of the C terms
                sum_over_zz = sum_over_zz + 2 * A(mm, dd) * mc(dd, nn) + ...
                                            2 * B(mm, dd) * ms(dd, nn);
            end
            normTerm = normTerm + y(mm, nn) ^ 2 ...
                                - 2 * y(mm, nn) * (sum_over_pp + C(mm)) ...
                                + sum_over_ij ...
                                + sum_over_qq ...
                                + C(mm) * (sum_over_zz + C(mm));
        end
    end


    % log p calculation
    log_p = N * (- sum(k + log(2 * pi * besseli(0, k, 1))) ...
                   + M / 2 * log(lambda2_y / (2 * pi)))...
            + prior - 0.5 * lambda2_y * normTerm;

    log_p = - log_p;

end

