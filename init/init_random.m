function [u, A, B, c, lambda2_y] = init_random(M, D)
%INIT_RADDOM Initialises the values of A, B, c, and lambda2_y using
%random numbers.
    A = rand(M, D);
    B = rand(M, D);
    u = rand(D, 1);
    c = rand(2, 1);
    lambda2_y = rand();
end