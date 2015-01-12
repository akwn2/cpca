function val = get_gvm_grid_moment1(k1, k2, m1, m2)
%GET_GVM_GRID_MOMENT0 Summary of this function goes here
%   Detailed explanation goes here
    fun = @(x, kappa1, kappa2, mu1, mu2) exp(1.i .* x) .* exp(kappa1 .* cos(x - mu1) + kappa2 .* cos(2 .* (x - mu2)));
    val = quadgk(@(x)fun(x, k1, k2, m1, m2), -pi, pi);
end

