function [mom0, mom1, mom2] = get_gvm_series_vectorized(mf_pars)
%GET_GVM_SERIES_VECTORIZED Summary of this function goes here
%   Detailed explanation goes here
    [k1, k2, m1, m2] = unpack_mf_pars(mf_pars);
    [n_dim, d_pts] = size(k1);
        
%     max_i = max(max(min(10 + ceil(0.5 * (k1 + k2)), 50)));
    max_i = 25;
    idx = 1:max_i;
    
    [~, ~, IDX] = meshgrid(zeros(1, d_pts), zeros(1, n_dim), idx);
    K1 = repmat(k1, [1 1 max_i]);
    K2 = repmat(k2, [1 1 max_i]);
    M1 = repmat(m1, [1 1 max_i]);
    M2 = repmat(m2, [1 1 max_i]);
    
    % Moment 0
    mom0 = besseli(0, k1, 1) .* besseli(0, k2, 1) ...
        + 2. * sum(besseli(2 .* IDX, K1, 1) .* cos(2. .* IDX .* (M1 - M2)) .* besseli(IDX, K2, 1), 3);
    
    mom0 = mom0 * 2. * pi;
    
    % Moment 1
    mom1 = besseli(1, k1, 1) .* besseli(0, k2, 1) ...
        + sum((exp(+2.0i .* IDX .* (M1 - M2)) .* besseli(2 .* IDX + 1, K1, 1) + ...
               exp(-2.0i .* IDX .* (M1 - M2)) .* besseli(2 .* IDX - 1, K1, 1)) .* besseli(IDX, K2, 1), 3);
           
    mom1 = mom1 .* 2. .* pi .* exp(1.0i .* m1);
    
    % Moment 2
    mom2 = besseli(2, k1, 1) .* besseli(0, k2, 1) ...
        + besseli(0, k1, 1) .* besseli(1, k2, 1) .* exp(-2.0i * (m1 - m2)) ...
        + sum(exp(+2.0i .* IDX .* (M1 - M2)) .* besseli(2 .* IDX + 2, K1, 1) .* besseli(IDX, K2, 1) + ...
              exp(-2.0j .* (IDX + 1) .* (M1 - M2)) .* besseli(2 .* IDX, K1, 1) .* besseli(IDX + 1, K2, 1), 3);
                    
    mom2 = mom2 .* 2. .* pi .* exp(2.0i .* m1);
    
    safety_tol = 0;
    idx = find(mom0 <= safety_tol);
    if ~isempty(idx)        
        % Try approximating moments of the GvM by von Mises when reasonable
        idx = find(mom0 <= safety_tol);
        if ~isempty(idx)
            ii = find(k1(idx) > 10 .* k2(idx));
            mom0(idx(ii)) = 2. * pi .* besseli(0, k1(idx(ii)), 1);
            mom1(idx(ii)) = 2. * pi .* besseli(1, k1(idx(ii)), 1) .* exp(1.i .* m1(idx(ii)));
            mom2(idx(ii)) = 2. * pi .* besseli(2, k1(idx(ii)), 1) .* exp(2.i .* m1(idx(ii)));
        end
        
        % Try approximating by a similar concentration ratio
        idx = find(mom0 <= safety_tol);
        if ~isempty(idx)
            max_idx = max(size(idx));
            order = 25;
            redux = 0.75;
            rk1 = redux * k1(idx);
            rk2 = redux * k2(idx);
            for ii = 1:max_idx
                mom0(idx(ii)) = get_gvm_series_moment0(rk1(ii), rk2(ii), m1(idx(ii)), m2(idx(ii)), order);
                mom1(idx(ii)) = get_gvm_series_moment1(rk1(ii), rk2(ii), m1(idx(ii)), m2(idx(ii)), order);
                mom2(idx(ii)) = get_gvm_series_moment2(rk1(ii), rk2(ii), m1(idx(ii)), m2(idx(ii)), order);
            end
        end
        
        idx = find(mom0 <= safety_tol);
        if ~isempty(idx)
            %If it is still wrong we'll try to calculate the expression
            %using grid integration as a last resort on scaled
            %concentrations
            ratio = k1(idx) ./ k2(idx);
            rk1 = 5 * ones(size(k1(idx)));
            rk2 = ratio .* rk1;
            max_idx = max(size(idx));
            for ii = 1:max_idx
                mom0(idx(ii)) = get_gvm_grid_moment0(rk1(ii), rk2(ii), m1(idx(ii)), m2(idx(ii)));
                mom1(idx(ii)) = get_gvm_grid_moment1(rk1(ii), rk2(ii), m1(idx(ii)), m2(idx(ii)));
                mom2(idx(ii)) = get_gvm_grid_moment2(rk1(ii), rk2(ii), m1(idx(ii)), m2(idx(ii)));
            end
        end
% 
%         % Try using a cruder approximation
%         idx = find(mom0 <= safety_tol);
%         if ~isempty(idx)
%             max_idx = max(size(idx));
%             for ii = 1:max_idx
%                 order = 25;
%                 while (mom0(idx(ii)) < 0)
%                     order = order - 1;
%                     mom0(idx(ii)) = get_gvm_series_moment0(k1(idx(ii)), k2(idx(ii)), m1(idx(ii)), m2(idx(ii)), order);
%                 end
%                 mom1(idx(ii)) = get_gvm_series_moment1(k1(idx(ii)), k2(idx(ii)), m1(idx(ii)), m2(idx(ii)), order);
%                 mom2(idx(ii)) = get_gvm_series_moment2(k1(idx(ii)), k2(idx(ii)), m1(idx(ii)), m2(idx(ii)), order);
%             end
%         end
        
        
        
        idx = find(mom0 <= safety_tol);
        if ~isempty(idx)
            % If we get to this point, then something really wrong happened
            % and we have to break and investigate it.
            disp('There is an exception in the moment calculations')
            keyboard;
            throw(MException('DO_VEM:ZeroNormalizingConstant', ...
                             'Normalizing constant evaluated to zero.'));
            
         end
    end
end
