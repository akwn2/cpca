function val = get_gvm_series_moment2(k1, k2, m1, m2, max_i)
%GET_GVM_GRID_ Calculate the 2nd moment of the GvM using series
%   Detailed explanation goes here
    if nargin < 5
       max_i = 25
    end
    idx = 1:max_i;
    val = besseli(2, k1, 1) * besseli(0, k2, 1) ...
        + besseli(0, k1, 1) * besseli(1, k2, 1) * exp(-2.0i * (m1 - m2)) ...
        + sum(exp(+2.0i .* idx .* (m1 - m2)) .* besseli(2 .* idx + 2, k1, 1) .* besseli(idx, k2, 1) + ...
                        exp(-2.0j .* (idx + 1) .* (m1 - m2)) .* besseli(2 .* idx, k1, 1) .* besseli(idx + 1, k2, 1));
    val = val * 2. * pi * exp(2.0i * m1);
end