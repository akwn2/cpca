%set_test_parameters
% Creates all test parameters for a CPCA problem:
%   N = number of data points
%   M = dimension of a single data point
%   D = number of hidden dimensions
%
%   Model: y = A * cos(\theta) + B * sin(\theta) + C + epsilon
%          epsilon = Normal(epsilon; 0, lambda2_y * eye(M));
%          theta = vonMises(kappa, mu);
%
%   Mean field distributions are GvM(k1, k2, m1, m2).
%   Expectation of trigonometric functions under the mean field (obtainable
%   by using the trigonometric moments of the GvM):
%       <cos> = mc, <sin> = ms, <cos * sin> = msc, <cos^2> = mc2, <sin^2> =
%       ms2;

clear; clc;

%% Set all parameters
%------------------------------
fprintf('Generating test parameters...');
tol = 1.E-6;
G = 1E4;
N = 1000;
M = 4;
D = 2;

theta = linspace(-pi, pi, G);
k1 = 1 .* ones(D, 1);
k2 = 0 .* ones(D, 1);
m1 = 0.5 .* ones(D, 1);
m2 = 0 .* ones(D, 1);

Z0 = zeros(D, 1);
Z1 = zeros(D, 1);
Z2 = zeros(D, 1);

mc = zeros(D, 1);
ms = zeros(D, 1);
mc2 = zeros(D, 1);
msc = zeros(D, 1);
ms2 = zeros(D, 1);

% Generate the "samples" from a Generalised von Mises
gvm = zeros(D, G);
T = zeros(D, N);
for dd = 1:D
    gvm(dd, :) = exp(k1(dd) .* cos(theta - m2(dd)) + ...
                     k2(dd) * cos(2 .* (theta - m2(dd))));
    T(dd,:) = randsample(theta, N, true, gvm(dd, :));
   
    % Zeroth moment (Partition function)
    f0 = @(x) exp(k1(dd) .* cos(x - m1(dd)) ...
                  + k2(dd) .* cos(2 .*(x - m2(dd))) ...
                  - k1(dd) - k2(dd));
    
    Z0(dd) = quadgk(f0, -pi, pi);
    
    % First moment
    f1 = @(x) exp(+1.j .* x + k1(dd) .* cos(x - m1(dd)) ...
                            + k2(dd) .* cos(2 .*(x - m2(dd))) ...
                            - k1(dd) - k2(dd));
    
    Z1(dd) = quadgk(f1, -pi, pi);
    
    % Second moment
    f2 = @(x) exp(+2.j .* x + k1(dd) .* cos(x - m1(dd)) ...
                            + k2(dd) .* cos(2 .*(x - m2(dd))) ...
                            - k1(dd) - k2(dd));

    Z2(dd) = quadgk(f2, -pi, pi);
end
            
M1 = repmat(Z1 ./ Z0 .* exp(1.j .* m1), [1, N]);
M2 = repmat(Z2 ./ Z0 .* exp(2.j .* m1), [1, N]);

mc = reshape(real(M1), D, N);
ms = reshape(imag(M1), D, N);
mc2 = reshape(0.5 .* (1 + real(M2)), D, N);
msc = reshape(0.5 .* imag(M2), D, N);
ms2 = reshape(0.5 .* (1 - real(M2)), D, N);

% Calculate the sine and cosine
cos_T = cos(T);
sin_T = sin(T);

A = [1, 0;...
     0, 0;...
     1, 1;...
     0, 0];
B = [0, 0;...
     1, 0;...
     0, 0;...
     1, 1];

c = [1.5, 2.5]';
C = [c', c']';
u = k1 .* cos(m1);
v = k1 .* sin(m1);
lambda2_y = 10.0;
rms = sqrt(1 / lambda2_y);

k1 = repmat(k1, [1, N]);
k2 = repmat(k2, [1, N]);
m1 = repmat(m1, [1, N]);
m2 = repmat(m2, [1, N]);

% Generate the data from the model
y = A * cos_T + B * sin_T + repmat(C, [1, N]) + rms .* randn(M, N);
fprintf(' Done.\n');
