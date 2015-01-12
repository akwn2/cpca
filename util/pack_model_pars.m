function model_pars = pack_model_pars(u, v, A, B, lambda2_y)
%PACK_MODEL_PARS Utility to pack model parameters
%   Detailed explanation goes here
    model_pars{1} = u;
    model_pars{2} = v;
    model_pars{3} = A;
    model_pars{4} = B;
    model_pars{5} = lambda2_y;
end

