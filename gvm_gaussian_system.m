function F = gvm_gaussian_system(varlist, k1, k2, m1, m2)
%GVM_SYSTEM System of equations to find the Gaussian corresponding to the
%GvM of the parameters supplied.
%   Detailed explanation goes here

    % unpack
    rho = varlist(1);
    sigma1 = varlist(2);
    sigma2 = varlist(3);
    nu1 = varlist(4);
    nu2 = varlist(5);

    % Compute constants
    c1 = k1 * cos(m1);
    c2 = k1 * sin(m1);
    c3 = k2 * cos(2 * m2);
    c4 = k2 * sin(2 * m2);
    c5 = k1^2;
    
    F = zeros(5,1);
%     F(1) = c1 - 1/(1 - rho^2) * (rho * nu2 / (sigma1 * sigma2) - nu1 / sigma1^2);
%     F(2) = c2 - 1/(1 - rho^2) * (rho * nu1 / (sigma1 * sigma2) - nu2 / sigma2^2);
%     F(3) = c3 - 1/(4 - 4 * rho^2) * (1 / sigma1^2 - 1 / sigma2^2);
%     F(4) = c4 - 1/(2 - 2 * rho^2) * (rho / (sigma1 * sigma2));
%     F(5) = c5 - (...
%         (1/(1 - rho^2) * (rho * nu2 / (sigma1 * sigma2) - nu1 / sigma1^2))^2 +...
%         (1/(1 - rho^2) * (rho * nu1 / (sigma1 * sigma2) - nu2 /
%         sigma2^2))^2);

    F(1) = c1 * (1 - rho^2) * (sigma1^2 * sigma2) - rho * nu2 * sigma1 - nu1 * sigma2;
    F(2) = c2 * (1 - rho^2) * (sigma1 * sigma2^2) - rho * nu1 * sigma2 - nu2 * sigma1;
    F(3) = c3 * (4 - 4 * rho^2) * (sigma1^2 * sigma2^2) - (sigma2^2 - sigma1^2);
    F(4) = c4 * (2 - 2 * rho^2) * (sigma1 * sigma2) - rho;
    F(5) = c5 - (F(1) ^ 2 + F(2) ^ 2);
end

