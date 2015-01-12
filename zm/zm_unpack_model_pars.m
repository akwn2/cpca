function [u, A, B, lambda2_y] = zm_unpack_model_pars(model_pars)
%UNPACK_MODEL_PARS Utility to unpack model parameters
%   Detailed explanation goes here
    u = model_pars{1};
    A = model_pars{2};
    B = model_pars{3};
    lambda2_y = model_pars{4};
end
