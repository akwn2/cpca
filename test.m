k1 = 100.0;
k2 = 100.0;

m1 = 0.0;
m2 = 0.0;

gvm_k = 1.0;
gvm_kpp = 0.0;

ii = 0;
while abs(gvm_kpp - gvm_k) > 1e-6
    ii = ii + 1;
    
    gvm_k = get_gvm_series_moment2(k1, k2, m1, m2, ii);
    gvm_kpp = get_gvm_series_moment2(k1, k2, m1, m2, ii + 1);
end
disp(ii)