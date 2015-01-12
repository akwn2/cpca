function val = get_gvm_series_moment1(k1, k2, m1, m2, max_i)
%GET_GVM_SERIES_MOMENTS1 Calculate the 1st moment of the GvM using series
%   Detailed explanation goes here
    if nargin < 5
       max_i = 25
    end
    
    idx = 1:max_i;
    val = besseli(1, k1, 1) * besseli(0, k2, 1) ...
        + sum((exp(+2.0i .* idx .* (m1 - m2)) .* besseli(2 .* idx + 1, k1, 1) + ...
               exp(-2.0i .* idx .* (m1 - m2)) .* besseli(2 .* idx - 1, k1, 1)) .* besseli(idx, k2, 1));
    val = val * 2. * pi * exp(1.0i * m1);
end
