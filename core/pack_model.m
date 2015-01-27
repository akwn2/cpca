function model = pack_model(u, v, A, B, lambda2_y)
%PACK_model Utility to pack model parameters
%   Detailed explanation goes here
    model{1} = u;
    model{2} = v;
    model{3} = A;
    model{4} = B;
    model{5} = lambda2_y;
end

