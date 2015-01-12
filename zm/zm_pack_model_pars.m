function model_pars = zm_pack_model_pars(u, A, B, lambda2_y)
%PACK_MODEL_PARS Utility to pack model parameters
%   Detailed explanation goes here
    model_pars{1} = u;
    model_pars{2} = A;
    model_pars{3} = B;
    model_pars{4} = lambda2_y;
end
