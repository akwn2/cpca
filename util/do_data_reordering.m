function [y, idx_map, comp_var] = do_data_reordering(y)
%DO_DATA_REORDERING Reorders data for algorithm to initalise and run
%   Detailed explanation goes here

    [m, ~] = size(y);
    aux = y;
    comps = m / 2;
    comp_var = zeros(comps, 1);
    
    % Calculate the variances of each pair of variables
    for ii = 1:comps
        comp_var(ii) = max(var(y(2 * ii - 1,:)), var(y(2 * ii,:)));
    end
    % Get the mapping from previous indexes to new indexes)
    [~, idx_map] = sort(comp_var);
    
    % Rearrange rows of y
    for ii = 1:comps
        y(2 * ii - 1,:) = aux(2 * idx_map(ii) - 1,:);
        y(2 * ii,:) = aux(2 * idx_map(ii),:);
    end
end

