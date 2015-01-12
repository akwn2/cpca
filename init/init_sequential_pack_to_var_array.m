function var_array = init_sequential_pack_to_var_array(u, v, A, B, lambda2_y)
%PACK_TO_VAR_ARRAY Packing to a single array for optimization
%   Detailed explanation goes here
    [dim_m, dim_n] = size(A);
    
    A = reshape(A, dim_m * dim_n, 1);
    B = reshape(B, dim_m * dim_n, 1);
    
    var_array = cat(1, u, v, A, B, lambda2_y);

end

