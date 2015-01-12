function mode = gvm_modes(k1, k2, m1, m2)
%GVM_MODES Finds the modes of the GvM(k1, k2, m1, m2) using its Gaussian
%construction.
%   Detailed explanation goes here
    f = @(x)gvm_gaussian_system(x,k1,k2,m1,m2);
    x = fsolve(f,rand(5,1));
    % Unpack
    rho = x(1);
    sigma1 = x(2);
    sigma2 = x(3);
    nu1 = x(4);
    nu2 = x(5);

    C = [sigma1^2, sigma1 * sigma2 * rho;
        sigma1 * sigma2 * rho, sigma2^2];

    % Line direction is the one given by the greatest eigenvector (max
    % variance => lowest decay rate)
    [V, D] = eig(C);
    d = diag(D);
    z = V(:,find(max(d)));

    a = z(1) / z(2);
    b = nu2 - nu1 * z(2) / z(1);
    
    delta = (a * b)^2 - 4 * (a ^ 2 + 1)*(b ^ 2 - 1);
    
    if delta < 0
        fprintf('negative delta!\n');
    else
        if abs(delta) < 1E-6
            mode = -b / (2 * a);
        else
            mode = zeros(1, 2);
            mode(1) = (-b + sqrt(delta)) / (2 * a);
            mode(2) = (-b - sqrt(delta)) / (2 * a);
        end
    end
    
end