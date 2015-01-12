function [mom0, mom1, mom2] = get_gvm_series_moments(mf_pars)
%GET_GVM_SERIES_MOMENTS Summary of this function goes here
%   Detailed explanation goes here
    [k1, k2, m1, m2] = unpack_mf_pars(mf_pars);
    
    [n_dim, d_pts] = size(k1);
    
    k1 = reshape(k1, n_dim * d_pts, 1);
    k2 = reshape(k2, n_dim * d_pts, 1);
    m1 = reshape(m1, n_dim * d_pts, 1);
    m2 = reshape(m2, n_dim * d_pts, 1);
    
    mom0 = arrayfun(@(x, y, z, w)get_gvm_series_moment0(x, y, x, w), k1, k2, m1, m2);
    mom1 = arrayfun(@(x, y, z, w)get_gvm_series_moment1(x, y, x, w), k1, k2, m1, m2);
    mom2 = arrayfun(@(x, y, z, w)get_gvm_series_moment2(x, y, x, w), k1, k2, m1, m2);
    
    mom0 = reshape(mom0, n_dim, d_pts);
    mom1 = reshape(mom1, n_dim, d_pts);
    mom2 = reshape(mom2, n_dim, d_pts);
end
