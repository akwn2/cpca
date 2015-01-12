function r = small_rand(n, m, w)
if nargin < 3
    w = 1E3;
end
    r = rand(n, m) / w;
end