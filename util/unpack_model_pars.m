function [u, v, A, B, lambda2_y] = unpack_model_pars(model_pars)
%UNPACK_MODEL_PARS Utility to unpack model parameters
%   Detailed explanation goes here
    u = model_pars{1};
    v = model_pars{2};
    A = model_pars{3};
    B = model_pars{4};
    lambda2_y = model_pars{5};
end

