function val = get_gvm_series_moment0(k1, k2, m1, m2, max_i)
%GET_GVM_SERIES_MOMENT0 Calculate the Zeroth moment of the GvM using series
%   Detailed explanation goes here
    if nargin < 5
       max_i = 25;
    end
    idx = 1:max_i;
    val = besseli(0, k1, 1) * besseli(0, k2, 1) ...
        + 2. * sum(besseli(2 .* idx, k1, 1) .* cos(2. .* idx .* (m1 - m2)) .* besseli(idx, k2, 1));
    val = val * 2. * pi;
end

