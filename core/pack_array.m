function vararray = pack_array(u, A, B, c, lambda2_y)
%PACK_ARRAY Packing to a single array for optimization
%   Detailed explanation goes here
    [M, D] = size(A);
    
    A = reshape(A, [M * D, 1]);
    B = reshape(B, [M * D, 1]);
    
    vararray = cat(1, u, A, B, c, lambda2_y);

end

