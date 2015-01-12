function [k1, k2, m1, m2] = unpack_mf_pars(mf_pars)
%UNPACK_MF_PARS Unpacks the mean field parameters
%   Detailed explanation goes here
    k1 = mf_pars{1};
    k2 = mf_pars{2};
    m1 = mf_pars{3};
    m2 = mf_pars{4};
end

