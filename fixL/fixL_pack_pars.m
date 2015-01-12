function pars = fixL_pack_pars(y, mc, ms, mc_sq, ms_sq, msc, n_dim, m_dim, d_pts, lambda2_y)
%PACK_PARS Utility to pack the averaged moments and data
%   Detailed explanation goes here
    pars{1} = y;
    pars{2} = mc;
    pars{3} = ms;
    pars{4} = mc_sq;
    pars{5} = ms_sq;
    pars{6} = msc;
    pars{7} = n_dim;
    pars{8} = m_dim;
    pars{9} = d_pts;
    pars{10} = lambda2_y;
end
