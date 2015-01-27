function [u, v, A, B, lambda2_y] = unpack_model(model)
%UNPACK_model Utility to unpack model parameters
%   Detailed explanation goes here
    u = model{1};
    v = model{2};
    A = model{3};
    B = model{4};
    lambda2_y = model{5};
end

