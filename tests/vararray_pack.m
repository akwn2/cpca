function vararray = vararray_pack(u, A, B, lambda2_y)
%VARARRAY_PACK Packing to a single array for optimization
%   Detailed explanation goes here
    [dim_m, dim_n] = size(A);
    
    A = reshape(A, dim_m * dim_n, 1);
    B = reshape(B, dim_m * dim_n, 1);
    
    vararray = cat(1, u, A, B, lambda2_y);

end

