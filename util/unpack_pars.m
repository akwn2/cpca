function [y, mc, ms, mc_sq, ms_sq, msc, n_dim, m_dim, d_pts] = unpack_pars(pars)
%UNPACK_PARS Utility to unpack the averaged moments and data
%   Detailed explanation goes here
    y = pars{1};
    mc = pars{2};
    ms = pars{3};
    mc_sq = pars{4};
    ms_sq = pars{5};
    msc = pars{6};
    n_dim = pars{7};
    m_dim = pars{8};
    d_pts = pars{9};
end

