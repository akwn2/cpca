function [k1, k2, m1, m2] = unpack_mf(mf)
%UNPACK_mf Unpacks the mean field parameters
%   Detailed explanation goes here
    k1 = mf{1};
    k2 = mf{2};
    m1 = mf{3};
    m2 = mf{4};
end

