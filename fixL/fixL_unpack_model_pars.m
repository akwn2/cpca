function [u, A, B] = fixL_unpack_model_pars(model_pars)
%UNPACK_MODEL_PARS Utility to unpack model parameters
%   Detailed explanation goes here
    u = model_pars{1};
    A = model_pars{2};
    B = model_pars{3};
end
